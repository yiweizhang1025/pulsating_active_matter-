/*

  This code simulates the dynamics of N active deforming particles:
  
  \dot r_i = - \mu_r \sum_{jvi} \nabla U(r_i-r_j, \phi_i, \phi_j) + \sqrt{2\mu_r T} \xi_i
	\dot\phi_i = \omega + \mu_\phi \sum_{jvi} [ \epsilon sin(\phi_i-\phi_j) - \partial_{\phi_i} U(r_i-r_j, \phi_i, \phi_j)] + \sqrt{2\mu_\phi T} \eta_i

	where

	range of interaction    = radius * [1 + lambda * sin(\phi)] / (1 + lambda)
	< \xi_i(t) \xi_j(0) >   = \delta_{ij} \delta(t)
	< \eta_i(t) \eta_j(0) > = \delta_{ij} \delta(t)


	Lx, Ly          = system size
	rho             = density
  T               = temperature
  mu_r            = mobility of position
  mu_phi          = mobility of size variable
  omega           = driving of size variable
  dt              = time step
  totaltime       = duration of the run
  initial_time    = waiting time before record
  interval_record = time interval between successive records
  interval_snap   = time interval between successive snapshots
  radius          = range of interaction
  amplitude_r     = strength of interaction
  amplitude_phi   = strength of alignment
	lambda          = interaction parameter
	dist            = radius of density averaging: 2*dist+radius < min(Lx,Ly)
	dx              = binning of density correlations
	drho            = binning of density
	dphi            = binning of size

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <time.h>
#include "mt19937-64.c"
#define PI 3.14159265358979323846

// Structure Particule: contains coordinates and box id of each particle
typedef struct {
  double qx; // Position
  double qy;
  double phi; // Size
  int i; // Box in which the particle is currently
  int j;
} particle;

// Structure Forces: contains x and y projection of interaction forces with neighboring particles
typedef struct {
  double fx;
  double fy;
  double torque;
} force;

// site = linked chained of all particle in a given box
// elem: itself + pointer to the next particles
typedef struct elem elem;

struct elem
{
  int num; // number of a particle in the box
  struct elem *next; // pointer onto the next element. next->num is the number of another particle in the box.
};
	
typedef elem* site; // Site is a pointer on a box, i.e. on a particle and on a pointer to another particle.


// Compute the sign of an argument
int sign(double x){
  return (x > 0) ? 1 : (-1);
}
int sign_int(int x){
  return (x > 0) ? 1 : (-1);
}

// Choose next and previous elements with periodic boundary condition
int NEXT(int k, double L){
  return (k == (int) L-1) ? 0 : (k+1);
}

int PREV(int k, double L){
  return (k == 0) ? ((int) L-1) : (k-1);
}

int IDX(int k, double L, int n){
  return (k == (int) L+n) ? 0 : (k+n);
}


/* ---------------------------
    DECLARATION OF FUNCTIONS
   ------------------------------*/

// Manipulate linked chain
site list_add(site list, int id); // Add a particle in a box
site list_remove(site list, int id); // Remove a particle from a box

// Store data
void Record(FILE* outputdata, particle *Particles, int N); // Particle coordinates
void Record_corr(FILE* output_corr, double dx, double *Correlation, long N_corr, double factor); // Density correlations
void Record_hist(FILE* output_hist, double dhist, double *Histogram, long N_hist, double factor); // Histogram array

// Start from particles randomly deposited in the system
void InitialConditions(int N, double Lx, double Ly, double radius, site **Sites, particle *Particles);

// Compute forces between neighboring particles
void add_force(int k1, particle p1, int k2, particle p2, double radius, double lambda, double amplitude_r, double amplitude_phi, double Ly, double Lx, force *Force);

// Implement the dynamics
void MoveParticles(int N, int Nx, int Ny, double Lx, double Ly, double *Noise_r, double *Noise_phi, double omega, double radius, double lambda, double amplitude_r, double amplitude_phi, double mu_r, double mu_phi, site **Sites, particle *Particles, double dt, double Time[0], force *Force);

// Analyze statistics
void Analysis_dens(int N, int Nx, int Ny, double Lx, double Ly, double radius, double lambda, site **Sites, particle *Particles, double **Density, int idx_max); // Local density
void add_density(particle p, double radius, double lambda, int i, int j, int idx_max, double Ly, double Lx, double **Density); // Compute overlap
void Analysis_corr(int N, double Lx, double Ly, double radius, particle *Particles, double dx, double *Correlation); // Density correlations


/* ---------------------------
    END DECLARATION OF FUNCTIONS
   ------------------------------*/

/* ---------------------------
    MAIN
   ------------------------------*/

int main(int argc, char *argv[]){ 
  /* 
     argc is the number of argument of the command line, and argv is a list
     of string containing the arguments. argv[0] is the name of the
     executable.
  */
  
  char params[] = "Input: output Lx Ly rho T mu_r mu_phi omega dt totaltime initial_time interval_record interval_snap seed radius amplitude_r amplitude_phi lambda dist drho dx dphi";

  // Check that the number of inputs is correct
  if(argc!=24){
    printf("%s\n",params);
    exit(1);
  }
  
  // DEFINITION OF VARIABLES  
  int i, j, k, idx; // Iterators
  int N; // Number of particles
  int Nx, Ny; // Number of boxes of length radius in Lx, Ly
  double Lx, Ly; // System size
  double rho; // Particle density
  double T; // Temperature
  double mu_r, mu_phi; // Mobilities
  double omega; // Driving of size
  double amplitude_r, amplitude_phi,amplitude_phi_max,amplitude_phi_min; // Strength of interactions
  double *Noise_r, *Noise_phi; // Array of noise terms
	double *Correlation; // Array of density correlations
	double **Density, *Histogram_dens; // Arrays of density
	double *Histogram_phi; // Array of size
	double drho, dx, dphi, rho_max, dist; // Binning of histogram
	long   N_corr, N_dens, N_phi; // Maximum value of histograms
	int    idx_max; // Size of density averaging in box length
	double mu_avg,mu2_avg;
	double Size, Size2, Current, dummy; // Order parameter and current of size
  double dt; // Time step
  double Time[0]; // Current time
  double totaltime; // Total duration of the run (after initialisation)
  double initial_time; // Time for equilibration
  double interval_record; // Time between two recordings
  double interval_snap; // Time between two snapshots
  double nb_record;  // Number of time steps of the next record
  double nb_snap;  // Number of time steps of the next snapshot
  long   counter_stat, aux; // Counter for statistics
  long   seed; // Seed of random number generator
  double radius,radius0, lambda; // Particle radius
  time_t clock; // Time to measure duration of the simulation
  
  FILE *output; // File where parameters are stored
  FILE *output_snap; // File where snapshots are stored
  FILE *output_corr; // File where correlation is stored
  FILE *output_dens; // File where density is stored
  FILE *output_phi; // File where size is stored
  FILE *output_order; // File where size order parameter is stored
  FILE *output_order2;
  char name[200]; // Name of output file containing the parameters
  char name_snap[200]; // Name of output file containing snapshots
  char name_corr[200]; // Name of output file containing correlation
  char name_dens[200]; // Name of output file containing density
  char name_phi[200]; // Name of output file containing size
  char name_order[200]; // Name of output file containing order parameter
  char name_order2[200];

  force *Force; // Forces applied on particles
  site **Sites; // Array of linked-chain corresponding to each boxes
  particle *Particles; // Particle coordinates and index

  // READING THE INPUT
  i = 1;
  sprintf(name, "%s", argv[i]);
  sprintf(name_snap, "%s_snap", name);
  sprintf(name_corr, "%s_corr", name);
  sprintf(name_dens, "%s_dens", name);
  sprintf(name_phi, "%s_phi", name);
  sprintf(name_order, "%s_order", name);
  sprintf(name_order2, "%s_order2", name);
  output      = fopen(name, "w");
  output_snap  = fopen(name_snap, "w");
  //output_corr  = fopen(name_corr, "w");
  //output_dens  = fopen(name_dens, "w");
  //output_phi   = fopen(name_phi, "w");
  output_order = fopen(name_order, "w");
  output_order2 = fopen(name_order2, "w");

  i++;
  Lx              = strtod(argv[i], NULL); i++;
  Ly              = strtod(argv[i], NULL); i++;
  rho             = strtod(argv[i], NULL); i++;
  T               = strtod(argv[i], NULL); i++;
  mu_r            = strtod(argv[i], NULL); i++;
  mu_phi          = strtod(argv[i], NULL); i++;
  omega           = strtod(argv[i], NULL); i++;
  dt              = strtod(argv[i], NULL); i++;
  totaltime       = strtod(argv[i], NULL); i++;
  initial_time    = strtod(argv[i], NULL); i++;
  interval_record = strtod(argv[i], NULL); i++;
  interval_snap   = strtod(argv[i], NULL); i++;
  seed            = strtol(argv[i], NULL, 10); i++;
  radius          = strtod(argv[i], NULL); i++;
  amplitude_r     = strtod(argv[i], NULL); i++;
  amplitude_phi_max   = strtod(argv[i], NULL); i++;
  amplitude_phi_min   = strtod(argv[i], NULL); i++;
  lambda          = strtod(argv[i], NULL); i++;
  dist            = strtod(argv[i], NULL); i++;
  drho            = strtod(argv[i], NULL); i++;
  dx              = strtod(argv[i], NULL); i++;
  dphi            = strtod(argv[i], NULL); i++;
	radius0=1.;

  // Print the parameters in output  
  fprintf(output, "%lg\n", Lx);
  fprintf(output, "%lg\n", Ly);
  fprintf(output, "%lg\n", rho);
  fprintf(output, "%lg\n", T);
  fprintf(output, "%lg\n", mu_r);
  fprintf(output, "%lg\n", mu_phi);
  fprintf(output, "%lg\n", omega);
  fprintf(output, "%lg\n", dt);
  fprintf(output, "%lg\n", totaltime);
  fprintf(output, "%lg\n", initial_time);
  fprintf(output, "%lg\n", interval_record);
  fprintf(output, "%lg\n", interval_snap);
  fprintf(output, "%lg\n", radius);
  fprintf(output, "%lg\n", amplitude_r);
  fprintf(output, "%lg\n", amplitude_phi_max);
  fprintf(output, "%lg\n", amplitude_phi_min);
  fprintf(output, "%lg\n", lambda);
  fprintf(output, "%lg\n", dist);
  fprintf(output, "%lg\n", drho);
  fprintf(output, "%lg\n", dx);
  fprintf(output, "%lg\n", dphi);
  fflush(output);

  // INITIALISATION OF VARIABLES
  init_genrand64(seed);
  
  // Number of particles given rho, Lx and Ly
  N = (int) floor(rho*Lx*Ly); 

  // Number of boxes of length radius in Lx, Ly
  Nx = (int) (Lx/radius0);
  Ny = (int) (Ly/radius0);  

  // Start the clock time
  clock = time(NULL);

  // Initially all linked chains are NULL pointers
  Sites = (site**) malloc(sizeof(site*)*Nx);
  for(i=0;i<Nx;i++){
    Sites[i] = malloc(sizeof(site)*Ny);
    for(j=0;j<Ny;j++)	Sites[i][j] = NULL; 
	}

  // Initialize arrays
  Density   = (double**) malloc(sizeof(double)*Nx);
  for(i=0;i<Nx;i++)	Density[i] = malloc(sizeof(double)*Ny);
  Noise_r   = (double*) malloc(sizeof(double)*2*N);
  Noise_phi = (double*) malloc(sizeof(double)*N);
  Particles = (particle*) malloc(sizeof(particle)*N);
  Force     = (force*) malloc(sizeof(force)*N);

  // Initialize density histogram
	rho_max        = 2*rho;
	N_dens         = (int) (rho_max/drho);
	idx_max        = (int) (dist/radius0);
  Histogram_dens = (double*) malloc(sizeof(double)*N_dens);
  memset(Histogram_dens, 0, N_dens*sizeof(double));

  // Initialize Density correlations
	if(8*radius0<Lx)	dist = 4*radius0;
	else 					  dist = Lx/2;
	N_corr      = (int) (dist/dx);
  Correlation = (double*) malloc(sizeof(double)*N_corr);
  memset(Correlation, 0, N_corr*sizeof(double));

  // Initialize size histogram
	N_phi         = (int) (2*PI/dphi);
  Histogram_phi = (double*) malloc(sizeof(double)*N_phi);
  memset(Histogram_phi, 0, N_phi*sizeof(double));

  // Initial conditions
  InitialConditions(N, Lx, Ly, radius0, Sites, Particles);
  printf("Random deposition succeeded\n");
  Time[0]      = -initial_time;
  Size         = 0;
  Current      = 0;
  nb_record    = 0;
  nb_snap      = 0;
  Size2		= 0;
  mu_avg	= 0;
  mu2_avg	= 0;
	counter_stat = 0;
	idx          = 0;
  output_order = fopen(name_order, "w");

   // Run dynamics until the time iterator reaches the final time
  while(Time[0]<totaltime){

    // Sample realizations of the noise 
    for(k=0;k<N;k++){
			Noise_r[2*k]   = sqrt(2*mu_r*T)*gasdevMT();
			Noise_r[2*k+1] = sqrt(2*mu_r*T)*gasdevMT();
			Noise_phi[k]   = sqrt(2*mu_phi*T)*gasdevMT();
		}
	if(Time[0]<0.) amplitude_phi = amplitude_phi_min; //fix the amplitude at min for relaxation
    if(Time[0]>0.) amplitude_phi = amplitude_phi_min + (amplitude_phi_max-amplitude_phi_min)*Time[0]/totaltime; //increment the amplitude from min
    // Move the particles according to the dynamics
	  MoveParticles(N, Nx, Ny, Lx, Ly, Noise_r, Noise_phi, omega, radius, lambda, amplitude_r, amplitude_phi, mu_r, mu_phi, Sites, Particles, dt, Time, Force);

		// Analysis of statistics
		if(idx==1){
			// Compute density
			//Analysis_dens(N, Nx, Ny, Lx, Ly, radius, lambda, Sites, Particles, Density, idx_max);

			// Compute density correlations						
			//Analysis_corr(N, Lx, Ly, dist, Particles, dx, Correlation);

			// Compute order parameter and current of size
			//dummy = 0;
		//	Size = 0;
		//	Current=0;
		//	for(i=0;i<N;i++){
		//		Current += omega + mu_phi*Force[i].torque;
		//		for(j=i+1;j<N;j++)	dummy += cos(Particles[i].phi-Particles[j].phi);
		//	}
		//	Size += sqrt(2*dummy);

			// Compute density histogram
			//for(i=0;i<Nx;i++){
			//	for(j=0;j<Ny;j++){
			//		aux = (int) (Density[i][j]/drho);
			//		Histogram_dens[aux]++;
				//}
			//}

			// Compute size histogram
		//	for(i=0;i<N;i++){
		//		aux = (int) ((Particles[i].phi - 2*PI*floor(Particles[i].phi/(2*PI)))/dphi);
		//		Histogram_phi[aux]++;
			//}

			// Increment of statistics counter
		//	counter_stat++;
		}

    // After thermalization
    if(Time[0]>=nb_record){
			//output_order = fopen(name_order, "w");
	 		// Print simulation progress
			printf("%s: %.4lg\n", name, 1e2*Time[0]/totaltime);

			// Allow record of statistics
			idx = 1;

			// Increase record counter
			nb_record += interval_record;

			if(Time[0]>=nb_snap){
				 if(idx==1){
                        // Compute density
                        //Analysis_dens(N, Nx, Ny, Lx, Ly, radius, lambda, Sites, Particles, Density, idx_max);

                        // Compute density correlations                                         
                        //Analysis_corr(N, Lx, Ly, dist, Particles, dx, Correlation);

                        // Compute order parameter and current of size
   	                	    	 dummy = 0;
        	               		 Size = 0;
					 Size2 = 0.;
                	       		// Current=0;
                      	 		 for(i=0;i<N;i++){
                               		//	 Current += omega + mu_phi*Force[i].torque;
                                		for(j=i+1;j<N;j++)      dummy += cos(Particles[i].phi-Particles[j].phi);
                    			    }
                       			 Size += sqrt(2*dummy);
					Size2 += 2*dummy;
                        // Compute density histogram
                        //for(i=0;i<Nx;i++){
                        //      for(j=0;j<Ny;j++){
                        //              aux = (int) (Density[i][j]/drho);
                        //              Histogram_dens[aux]++;
                                //}
                        //}

                        // Compute size histogram
                //      for(i=0;i<N;i++){
                //              aux = (int) ((Particles[i].phi - 2*PI*floor(Particles[i].phi/(2*PI)))/dphi);
                //              Histogram_phi[aux]++;
                        //}

                        // Increment of statistics counter
                   	     counter_stat++;
               		 }

				// Record statistics
				//remove(name_corr);
				//remove(name_dens);
				//remove(name_phi);
				//remove(name_order);

				//output_corr  = fopen(name_corr, "w");
				//output_dens  = fopen(name_dens, "w");
				//output_phi   = fopen(name_phi, "w");
				//output_order = fopen(name_order, "w");

				//Record_corr(output_corr, dx, Correlation, N_corr, PI*N*rho*counter_stat/4);
				//Record_hist(output_dens, drho, Histogram_dens, N_dens, counter_stat);
				//Record_hist(output_phi, dphi, Histogram_phi, N_phi, N*counter_stat);

				fprintf(output_order, "%.15lg\t%.15lg\n",amplitude_phi, 1-Size/N);
				fflush(output_order);
				fprintf(output_order2, "%.15lg\t%.15lg\n",amplitude_phi, Size2/(N*N));
                                fflush(output_order2);
				//fclose(output_order);
				mu_avg += Size/N;
				mu2_avg += Size2/(N*N);
				// Record snapshot of system
			//	Record(output_snap, Particles, N);

				// Increase snap counter
				nb_snap += interval_snap;
				
			}
		}
  }

  // Return the duration of the simulation
  printf("Simulation duration: %ld seconds\n", (long) time(NULL)-clock);
  fprintf(output, "%ld\n", (long) time(NULL)-clock);
  fprintf(output, "%ld\n", counter_stat);
  fprintf(output, "%f\n", mu_avg);
  fprintf(output, "%f\n", mu2_avg);
  fclose(output_order);
  return 0;
}

/* ---------------------------
    END MAIN
   ------------------------------*/

/* ---------------------------
    DYNAMICS
   ------------------------------*/

void MoveParticles(int N, int Nx, int Ny, double Lx, double Ly, double *Noise_r, double *Noise_phi, double omega, double radius, double lambda, double amplitude_r, double amplitude_phi, double mu_r, double mu_phi, site **Sites, particle *Particles, double dt, double Time[0], force *Force){
  
  int k1, k2, k, i, j, i_new, j_new;
  double dqx, dqy, dphi; // Increment of coordinates
	double test, max_amp=0; // Adaptative time stepping
	double dt_g, sdt_g; // Time step
	double radius0=1.;
  particle p1, p2;
  site site1; // Iterator on the particles of the site (i,j)
  site site2; // Iterator on the particles of the neighboring sites of the site (i,j)

  // Reset interaction forces
  memset(Force, 0, N*sizeof(force));
	
	// Iterate over neighboring particles
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      site1 = Sites[i][j]; // Pointer on particles in site (i,j)
      
      while(site1!=NULL){ // Iteration on particles in site (i,j)

				k1 = site1->num; // Number of particle 1
				p1 = Particles[k1]; // Structure of particle 1
	
				// Pointer on the next particle in the SAME site (i,j) 
				site2 = site1->next;
				while(site2!=NULL){ // While there is such a particle
					k2 = site2->num; // Get its number
					p2 = Particles[k2]; // Get its structure

					// Compute forces from particles on the SAME site
					add_force(k1, p1, k2, p2, radius, lambda, amplitude_r, amplitude_phi, Ly, Lx, Force); 
					site2 = site2->next;
				}

				// Site on the right
				site2 = Sites[NEXT(i,Nx)][j];
				while(site2!=NULL){
					k2 = site2->num;
					p2 = Particles[k2];

					// Compute forces from neighboring particles
					add_force(k1, p1, k2, p2, radius, lambda, amplitude_r, amplitude_phi, Ly, Lx, Force); 
					site2 = site2->next;
				}
	
				// Site on the bottom left
				site2 = Sites[PREV(i,Nx)][PREV(j,Ny)];
				while(site2!=NULL){
					k2 = site2->num;
					p2 = Particles[k2];

					// Compute forces from neighboring particles
					add_force(k1, p1, k2, p2, radius, lambda, amplitude_r, amplitude_phi, Ly, Lx, Force); 
					site2 = site2->next;
				}
	
				// Site on the bottom 
				site2 = Sites[i][PREV(j,Ny)];
				while(site2!=NULL){
					k2 = site2->num;
					p2 = Particles[k2];

					// Compute forces from neighboring particles
					add_force(k1, p1, k2, p2, radius, lambda, amplitude_r, amplitude_phi, Ly, Lx, Force); 
					site2 = site2->next;
				}
	
				// Site on the bottom right
				site2 = Sites[NEXT(i,Nx)][PREV(j,Ny)];
				while(site2!=NULL){
					k2 = site2->num;
					p2 = Particles[k2];

					// Compute forces from neighboring particles
					add_force(k1, p1, k2, p2, radius, lambda, amplitude_r, amplitude_phi, Ly, Lx, Force); 
					site2 = site2->next;
				}

				site1 = site1->next;
      }
    }
  }

	// Compute maximum of interaction amplitude
	for(k=0;k<N;k++){
		test = mu_r*fabs(Force[k].fx)/radius;
		if(test>max_amp) max_amp = test;
		test = mu_r*fabs(Force[k].fy)/radius;
		if(test>max_amp) max_amp = test;
		test = (omega + mu_phi*Force[k].torque)/(2*PI);
		if(test>max_amp) max_amp = test;
	}

 	// Adaptative time stepping
	if(max_amp*dt<1e-1)	dt_g = dt;
	else								dt_g = 1e-1/max_amp;
	sdt_g    = sqrt(dt_g);
	Time[0] += dt_g;	

  // Update coordinates of all particles
  for(k=0;k<N;k++){ // Iterations over the N particles
    
    // Extract position and label of the box containing the particle
    i = (int) (Particles[k].qx/radius0);
    j = (int) (Particles[k].qy/radius0);

	  // Dynamics
		dqx  = dt_g*mu_r*Force[k].fx + sdt_g*Noise_r[2*k];
	  dqy  = dt_g*mu_r*Force[k].fy + sdt_g*Noise_r[2*k+1];
		dphi = dt_g*(omega + mu_phi*Force[k].torque) + sdt_g*Noise_phi[k];

    // New coordinates after one iteration
    Particles[k].qx  += dqx;
    Particles[k].qy  += dqy;
    Particles[k].phi += dphi;

    // Pay attention to the boundary conditions
    if(Particles[k].qx>Lx || Particles[k].qx<0)  Particles[k].qx -= Lx*floor(Particles[k].qx/Lx);
		if(Particles[k].qy>Ly || Particles[k].qy<0)  Particles[k].qy -= Ly*floor(Particles[k].qy/Ly);

    // New label of the box containing the particle
    i_new = (int) (Particles[k].qx/radius0);
    j_new = (int) (Particles[k].qy/radius0);

    // Update de Sites
    if(i!=i_new || j!=j_new){
      Sites[i][j]         = list_remove(Sites[i][j], k);
      Sites[i_new][j_new] = list_add(Sites[i_new][j_new], k);
    }
  }
}

/* ---------------------------
    END DYNAMICS
   ------------------------------*/

/* ---------------------------
    DATA ANALYSIS
   ------------------------------*/

void Analysis_corr(int N, double Lx, double Ly, double dist, particle *Particles, double dx, double *Correlation){

	int i, j;
	long aux; // Auxiliary variable
	double r, rx, ry; // Interparticle distance

	for(i=0;i<N;i++){
		for(j=i+1;j<N;j++){
			// Interparticle distance
			rx = Particles[i].qx - Particles[j].qx;
			ry = Particles[i].qy - Particles[j].qy;

			// Pay attention to boundary conditions
			rx -= Lx*floor(rx/Lx);
			ry -= Ly*floor(ry/Ly);
			r   = sqrt(rx*rx + ry*ry);

			// Check if position is in circle of radius 'dist'
			if(r<dist){
				// Position in histogram array
				aux = (int) (r/dx);
				Correlation[aux]++;
			}
		}
	}
}

void Analysis_dens(int N, int Nx, int Ny, double Lx, double Ly, double radius, double lambda, site **Sites, particle *Particles, double **Density, int idx_max){

	int i, j, ii, jj, idx, jdx;
  particle p; // Iterator on particles in sites
  site site; // Iterator on sites in grid

	for(i=0;i<Nx;i++){
	  for(j=0;j<Ny;j++){
			// Initialize density
			Density[i][j] = 0;

			for(ii=0;ii<2*idx_max+1;ii++){
				for(jj=0;jj<2*idx_max+1;jj++){
					// Compute site indices with boundary conditions
					idx = i+ii-idx_max;
					jdx = j+jj-idx_max;
					if(idx<0 || idx>Nx-1)	idx -= Nx*sign_int(idx);
					if(jdx<0 || jdx>Ny-1)	jdx -= Ny*sign_int(jdx);

					// Pointer on particles in site (idx,jdx)
			    site = Sites[idx][jdx];

			    while(site!=NULL){
						// Select next particle in the box
						p = Particles[site->num];

						// Compute overlap between particles
						add_density(p, radius, lambda, i, j, idx_max, Ly, Lx, Density);
						site = site->next;
					}
				}
			}
		}
	}
}

// Compute the overlap between particles and add contribution to density array
void add_density(particle p, double radius, double lambda, int i, int j, int idx_max, double Ly, double Lx, double **Density){

	double rad, rad2, dist, dist2, rx, ry, r, r2; // Distance from box center
	double area = 0; // Overlap area

	// Distance from box center
	rx   = (i+.5)*radius - p.qx;
	ry   = (j+.5)*radius - p.qy;
	rad  = radius/2/(1+lambda)*(1 + lambda*sin(p.phi));
	dist = (idx_max+.5)*radius;

  // Pay attention to boundary conditions
  if(fabs(rx)>(Lx-rad-dist))	rx -= sign(rx)*Lx;
  if(fabs(ry)>(Ly-rad-dist))	ry -= sign(ry)*Ly;
	r = sqrt(rx*rx + ry*ry);

	// Compute overlap
	if(r<dist+rad){
		if(r<dist-rad)	area = PI*rad*rad; // Complete overlap
		else{
			// Compute distances
			r2    = r*r;
			rad2  = rad*rad;
			dist2 = dist*dist;

			// Partial overlap
			// Ref: https://mathworld.wolfram.com/Circle-CircleIntersection.html
			area  = rad2*acos((r2+rad2-dist2)/(2*r*rad));
			area += dist2*acos((r2-rad2+dist2)/(2*r*dist));
			area -= .5*sqrt((-r+rad-dist)*(-r-rad+dist)*(-r+rad+dist)*(r+rad+dist));
		}
	}

	// Add contribution to density
	Density[i][j] += area/(PI*dist*dist);
}


/* ---------------------------
    END DATA ANALYSIS
   ------------------------------*/

/* ---------------------------
    FORCES
   ------------------------------*/

// Compute the force between particles and add it to the force array
void add_force(int k1, particle p1, int k2, particle p2, double radius, double lambda, double amplitude_r, double amplitude_phi, double Ly, double Lx, force *Force){

  double rad, rad6, rx, ry, r, r6; // Interparticle distance
  double amp_rep, amp_phi, amp_align; // Interaction strengths
  
  // Interparticle distance
  rx  = p2.qx - p1.qx;
  ry  = p2.qy - p1.qy;
	rad = radius/(1+lambda)*(1 + lambda/2*(sin(p1.phi) + sin(p2.phi)));

  // Pay attention to boundary conditions
  if(fabs(rx)>(Lx-rad)) rx -= sign(rx)*Lx;
  if(fabs(ry)>(Ly-rad)) ry -= sign(ry)*Ly;
	r = sqrt(rx*rx + ry*ry);

	// Add contribution to the force
	if(r<rad){
		// Interaction strengths
		r6        = r*r*r*r*r*r;
		rad6      = rad*rad*rad*rad*rad*rad;
		amp_rep   = 12*amplitude_r*rad6/r6/r*(rad6/r6 - 1);
//		amp_rep   = 2*amplitude_r*(1 - r/rad)/rad;
		amp_phi   = amp_rep*radius*r*lambda/(2*rad*(1+lambda));
		amp_align = amplitude_phi*sin(p2.phi-p1.phi);		

		// Interaction on p1
		Force[k1].fx     -= amp_rep*rx/r;
		Force[k1].fy     -= amp_rep*ry/r;
		Force[k1].torque += amp_align - amp_phi*cos(p1.phi);
		
		// Interaction on p2
		Force[k2].fx     += amp_rep*rx/r;
		Force[k2].fy     += amp_rep*ry/r;
		Force[k2].torque -= amp_align + amp_phi*cos(p2.phi);
	}
}

/* ---------------------------
    END FORCES
   ------------------------------*/

/* ---------------------------
    INITIAL CONDITIONS
   ------------------------------*/

// Deposition of N particles uniformly in the system
void InitialConditions(int N, double Lx, double Ly, double radius, site **Sites, particle *Particles){

  int k = 0; // Iterator on particles
  int i, j; // Indices of the box
  double X, Y; // Sampled coordinates

  // Iterations over the deposited particles
  while(k<N){
    
		// Initialization
		X = genrand64_real2()*Lx;
		Y = genrand64_real2()*Ly;

    // Extract label (i,j) of corresponding box
    i = (int) (X/radius);
    j = (int) (Y/radius);
  
    // Store spatial coordinates
    Particles[k].qx = X;
    Particles[k].qy = Y;

    // Store lattice coordinates
    Particles[k].i = i;
    Particles[k].j = j;
  
    // Initial size variables
    Particles[k].phi = genrand64_real2()*2*PI;
 
    // Add particle in the system
    Sites[i][j] = list_add(Sites[i][j], k);
    k++;
  }
}

/* ---------------------------
    END INITIAL CONDITIONS
   ------------------------------*/

/* ---------------------------
    MANIPULATE LINKED CHAINS
   ------------------------------*/

// Add new element at the beginning of the linked chain
site list_add(site list, int num){
  site new_elem  = malloc(sizeof(elem)); // Define temporary element

  new_elem->num  = num; // Assign label *num* to the temporary element
  new_elem->next = list; // Assign pointer of the temporary element to the linked chain *list* 

  return new_elem;
}

// Remove element k from linked chain
site list_remove(site list, int k){
  site it, it_prec;
  it = list;

  if(it->num==k){ // If element k is the first element
    free(list); // Remove the first element
    return it->next;
  }
 
  else{
    it_prec = list;
    it      = it->next;

    while(it!=NULL){ // While there is an element in the list

      if(it->num==k){
				it_prec->next = it->next;
				free(it); // Remove the element
				return list;
      }

      it_prec = it;
      it      = it->next;
    }
  }
  printf("Error: element is not found\n");
}

/* ---------------------------
    END MANIPULATE INKED CHAINS
   ------------------------------*/

/* ---------------------------
    PRINT
   ------------------------------*/

// Record particle coordinates
void Record(FILE* outputdata, particle* Particles, int N){

  int i;
  
  for(i=0;i<N;i++){
    fprintf(outputdata, "%lg\t%lg\t%lg\n", Particles[i].qx, Particles[i].qy, Particles[i].phi);
  }

  fflush(outputdata);
}

// Record correlation
void Record_corr(FILE* output_corr, double dx, double *Correlation, long N_corr, double factor){

  int i; // Iterator
  
	// Print data
  for(i=0;i<N_corr;i++) fprintf(output_corr, "%lg\t%lg\n", dx*(i+.5), Correlation[i]/(dx*dx*(i+.5)*factor));

  fflush(output_corr);
  fclose(output_corr);
}

// Record histogram
void Record_hist(FILE* output_hist, double dhist, double *Histogram, long N_hist, double factor){

  int i; // Iterator
  
	// Print data
  for(i=0;i<N_hist;i++) fprintf(output_hist, "%lg\t%lg\n", dhist*(i+.5), Histogram[i]/(dhist*factor));

  fflush(output_hist);
  fclose(output_hist);
}
/* ---------------------------
    END PRINT
   ------------------------------*/


