#ifndef ODESOLVER_H
#define ODESOLVER_H       //after this is defined the first time, compiler will
                          //skip over header if its included again b/c #ifndef (if not defined) will fail

#include "planet.h"
#include <vector>
#include <fstream>

using std::vector;

class ODEsolver
{
public:

    // properties
    double Gconst;
    int total_planets;
    vector<planet> all_planets;   //vector with  elements of objects planets
    double totalKinetic;
    double totalPotential;

    //Constructor
    ODEsolver(int zero);

    // functions
    void add(planet newplanet);
    void print_position(std::ofstream &output, double time, int number);
    void print_energy(std::ofstream &output, double time);
    void Euler(int IntegrationPoints, double FinalTime);
    void VelocityVerlet(int dimension, int integration_points, double final_time, int print_number, double epsilon);
    void GravitationalForce(planet &current, planet &other, double &Fx, double &Fy, double &Fz, double epsilon);
    void KineticEnergySystem();
    void PotentialEnergySystem();
    double EnergyLoss();
    bool Bound(planet OnePlanet);
};

#endif // ODESOLVER_H
