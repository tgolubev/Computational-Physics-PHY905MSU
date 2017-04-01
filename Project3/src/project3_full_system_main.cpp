//This is for the full solar system,  including all planets.

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>  //for files
#include <random>
#include <chrono>
#include <time.h>
//#include "planet.h"   odesolver.h already includes planet.h
#include "odesolver.h"

using namespace std;
using namespace chrono;

int main()
{
    // Numerical setup
    int integration_points = 50000;  // No. of integration points
    double final_time = 300.;       // End time of calculation (150  is not long enough to complete 1 orbit of all planets)
    bool corrections = false;       //No relativistic corrections applied

    //Set-up planets, Using April 1st, 2017 positions and velocities (in Au and Au/year), use 1yr = 365.25days
    //These numbers are based on origin of coordinate system at center of mass of solar system.
    //i.e. planet1 (name,mass,x,y,z,vx,vy,vz), name must be in " " marks
    planet planet1("Sun",1., 3.109182869400257e-03, 4.536032581431210e-03, -1.502332628558876e-04, -0.00130075, 0.00239436, 0.000028737);
    planet planet2("Mercury",1.66012e-7,-1.997697031204302e-01, 2.596435954332266e-01, 3.930764919716435e-02,-10.1211, -5.98857, 0.438896);
    planet planet3("Venus",2.44781e-6,-6.914113410296975e-01,-1.859648579125333e-01, 3.731440216587294e-02,1.9042, -7.15541, -0.208101);
    planet planet4("Earth",3.003467e-6,-9.770327933629532E-01, -1.899093075015161E-01, -1.346472796698014E-04,1.12159, -6.18619, 0.000269565);
    planet planet5("Mars",3.22713e-7,5.871508530860990e-01, 1.400379469087186, 1.476688515437949e-02,-4.52266, 2.41061, 0.161458);
    planet planet6("Jupiter",9.5458e-4,-5.199346663609341,-1.634983952396202, 1.230695166574554e-01,0.794561, -2.49878, -0.00740304);
    planet planet7("Saturn",2.85812e-4,-1.397265846418547, -9.948925301628403, 2.285941805683210e-01, 1.90617, -0.289666, -0.0709435);
    planet planet8("Uranus",4.36576e-5,1.819829127353053e+01,  8.138008397092724, -2.055373112007303e-01,-0.596947, 1.24444, 0.0123422);
    planet planet9("Neptune",5.1503e-5,2.842772334486519e+01, -9.420046539130526, -4.611579863476090e-01,0.352966, 1.09503, -0.0307882);
    planet planet10("Pluto",6.583e-9,9.939352708227002, -3.177412405108472e+01,  5.249692534427017e-01,1.11622, 0.0998563, -0.337162);

    //Setup the solar system system
    ODEsolver solar_system;         //create object of class ODEsolver with default constructor ODEsolver(). If put the () in declaration here, it doesn't work!
    solar_system.add(planet1);      //add planets to the solver
    solar_system.add(planet2);
    solar_system.add(planet3);
    solar_system.add(planet4);
    solar_system.add(planet5);
    solar_system.add(planet6);
    solar_system.add(planet7);
    solar_system.add(planet8);
    solar_system.add(planet9);
    solar_system.add(planet10);


    //Output the properties of the planets
    for(int i=0;i<solar_system.total_planets;i++)
    {
      cout << solar_system.all_planets[i].name << "'s Mass = " <<solar_system.all_planets[i].mass<< endl;
      cout << solar_system.all_planets[i].name <<"'s Initial Position = " << solar_system.all_planets[i].position[0] << ","
           << solar_system.all_planets[i].position[1]<<","<< solar_system.all_planets[i].position[2]<< endl;
    }
    //Tests of the setup
    cout << "Gconst = " <<solar_system.Gconst << endl;
    cout << "Number of Planets = " <<solar_system.total_planets<<endl;
    /*test planet_names (vector of strings)
    //for(int i=0; i<binary.total_planets;i++){
    //cout << "Planet" << i+1 << "'s name is " << binary.planet_names[i] <<endl;
    }*/

    /*
    //Euler method
    //start clock timer
    high_resolution_clock::time_point start1 = high_resolution_clock::now();

    solar_system.Euler(integration_points, final_time, corrections);  //Run Euler's method ODEsolver

    //stop clock timer and output time duration
    high_resolution_clock::time_point finish1 = high_resolution_clock::now();
    duration<double> time1 = duration_cast<duration<double>>(finish1-start1);
    cout << "Euler Solver CPU time = " << time1.count() << endl;

    */

    // Velocity Verlet
    //start clock timer
    high_resolution_clock::time_point start2 = high_resolution_clock::now();

    solar_system.VelocityVerlet(integration_points, final_time, corrections); //Run VVerlet ODEsolver

    //stop clock timer and output time duration
    high_resolution_clock::time_point finish2 = high_resolution_clock::now();
    duration<double> time2 = duration_cast<duration<double>>(finish2-start2);
    cout << "Velocity Verlet Solver CPU time = " << time2.count() << endl;


    /*  // RK4
        solver binary_rk(5.0);
        binary_rk.add(planet1);
        binary_rk.add(planet2);

        for(int j=0;j<dimension;j++){
             x[j] = planet1.position[j];
             v[j] = planet1.velocity[j];
         }
     */

    return 0;
}

