class Particle{
private:
    // Private members can be added here if needed
    int id; // Particle ID
    double r;
    // Particle motion
    double x, y; // Particle coordinates
    double vx, vy; // Particle velocities

    double force; // Force acting on the particle, if needed

public:
    // Constructor
    void Particle(int id, double r, double x, double y, double vx = 0.0, double vy = 0.0){
        this->id = id;
        this->r = r; 
        this->x = x;
        this->y = y;
        this->vx = vx;
        this->vy = vy;
    };
    
    // Settors
    void SetCoordinates(double x, double y){
        this->x = x;
        this->y = y;
    };

    void SetVelocity(double vx, double vy){
        this->vx = vx;
        this->vy = vy;
    }:

    void SetRadius(double r){
        this->r = r;
    };

    void SetForce(double force){
        this->force = force;
    };

    // Gettors
    int GetID() const {
        return id;
    };
    
    double GetRadius() const {
        return r;
    };

    double GetX() const {
        return x;
    };

    double GetY() const {
        return y;
    }; 
    
    double GetVx() const {
        return vx;
    };

    double GetVy() const {
        return vy;
    };  

    double GetForce() const {
        return force;
    };
}