#include <vector>

class Site {
private:
    int idX, idY; // Coordinates of the site in the grid
    vector<int> particleIDs; 
        vector<int> getParticleIDs(){
        return particleIDs;
    }
public:
    void popParticleID(int id){
        for (int k = 0; k < this->particleIDs.size(); ++k){
            if(particleIDs[k]==id){
                particleIDs.erase(particleIDs.begin() + k);
                break; // Exit after removing the first occurrence
            }
        }
        throw runtime_error("Error: particle ID not found in site");
    }
    /*
    void copy_particle_ids_to_deque(deque<int> &d){
        for(auto id : this->particle_ids){
        d.push_back(id);
        }
    }
    */
    
    void addParticleID(int id){
        for (int k = 0; k < this->particleIDs.size(); ++k){
            if(particleIDs[k]==id){
                throw runtime_error("Error: particle ID already exists in site");
            }
        }
        this->particleIDs.push_back(id);
    }
};