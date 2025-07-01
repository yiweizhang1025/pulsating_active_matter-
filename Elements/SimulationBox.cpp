
#include <string>

class SimulationBox {
    private:
        // Private member variables
        int width;
        int height;
        
        int nbOfParticles; // Number of particles in the simulation box

        double gridSpacing; // Spacing between grid points

        vector<vector<Site>> sites; // Density distribution in the simulation box
        vector<Particle> particles; // Particles

        string initializationType = "Random"; // Type of initialization (e.g., random, grid)
        string boundaryCondition = "Periodic"; // Boundary condition type
    public: 
        // Constructor
        SimulationBox(int width, int height, int nbOfParticles, double gridSpacing, string& initializationType = "Random", string& boundaryCondition = "Periodic") {
            this->width = width;
            this->height = height;
            this->nbOfParticles = nbOfParticles;
            this->gridSpacing = gridSpacing;
            this->initializationType = initializationType;
            this->boundaryCondition = boundaryCondition;
            // Initialize the sites grid
            int gridWidth = static_cast<int>(width / gridSpacing);
            int gridHeight = static_cast<int>(height / gridSpacing);
            sites.resize(gridWidth, vector<Site>(gridHeight));  
            particles.resize(nbOfParticles); 
        }

        // Getters
        int getWidth() const { return width; }
        int getHeight() const { return height; }
        int getNbOfParticles() const { return nbOfParticles; }
        double getGridSpacing() const { return gridSpacing; }
        vector<vector<Site>> getSites() const { return sites; }
        string getInitializationType() const { return initializationType; }
        string getBoundaryCondition() const { return boundaryCondition; }

        // Setters
        void setInitializationType(const string& type) { initializationType = type; }
        void setBoundaryCondition(const string& condition) { boundaryCondition = condition; }

        // Spread particles in the simulation box
        void spreadParticles() {
            // Implementation of particle spreading based on initializationType
            if (initializationType == "Random") {
                // Randomly distribute particles in the simulation box
                for (int i = 0; i < nbOfParticles; ++i) {
                    double x = static_cast<double>(rand()) / RAND_MAX * width;
                    double y = static_cast<double>(rand()) / RAND_MAX * height;
                    particles[i] = Particle(i, 1.0, x, y); // Assuming radius of 1.0 for simplicity
                    int gridX = static_cast<int>(x / gridSpacing);
                    int gridY = static_cast<int>(y / gridSpacing);
                    sites[gridX][gridY].addParticleID(i);
                }
            } else {
                throw runtime_error("Unknown initialization type: " + initializationType);
            }
        };
};