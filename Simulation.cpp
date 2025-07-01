#include "Elements/SimulationBox.cpp"

class Simulation{
    private:
        double Time; // Current simulation time
        double timeStep;

        double nbOfParticles; // Number of particles in the simulation
        double boxSizeX, boxSizeY; // Size of the simulation box
        double particleRadius; // Radius of the particles

        SimulationBox box; // Simulation box object

        string forcePotential;
        // Add more members as needed
        string boundaryCondition; // Type of boundary condition (e.g., periodic, reflective)
        string initialCondition; // Initial condition of the particles
        string outputFileName; // Name of the output file
    public:
        Simulation(double Time, double timeStep, double boxSizeX, double boxSizeY, double particleRadius, double nbOfParticles, string outputFileName = "simulation_output.txt", string boundaryCondition = "Periodic", string initialCondition = "Random", string forcePotential = "Lennard-Jones"){
            this ->Time = Time;
            this ->timeStep = timeStep;
            this ->boxSizeX = boxSizeX;
            this ->boxSizeY = boxSizeY;
            this ->particleRadius = particleRadius;
            this ->nbOfParticles = nbOfParticles;
            this ->boundaryCondition = "periodic"; // Default boundary condition
            this ->initialCondition = "random"; // Default initial condition    
            this ->outputFileName = outputFileName; // Default output file name
            this ->forcePotential = forcePotential; // Default force potential
        }
        
        void initializeBox(){
            box = SimulationBox(boxSizeX, boxSizeY, particleRadius, boundaryCondition);
            box.spreadParticles(); // Spread particles in the box based on initial condition
        }

        void runSimulation(){
            // Main simulation loop
            for (double t = 0; t < Time; t += timeStep) {
                // Update particle positions and velocities based on forces
                // Apply boundary conditions
                // Output results to file
            }
        }
}