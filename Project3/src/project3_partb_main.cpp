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
    int integration_points = 10000;  // No. of integration points
    double final_time = 50.;       // End time of calculation

    //Set-up planets
    planet planet1("Sun",1.,0.,0.,0.,0.,0.,0.);              // planet1 (name,mass,x,y,z,vx,vy,vz)  //NOTE: right now the Sun is at origin of coordinate system
                                                             //later will set the true COM of solar system to be the origin.
    planet planet2("Earth",0.000003,1.,0.0,0.0,0.0,6.3,0.);  // planet2 (name,mass,x,y,z,vx,vy,vz)   //name must be in " " marks

    //Output the properties of the planets
    cout << planet1.name << "'s Mass = " <<planet1.mass<< endl;
    cout << planet1.name <<"'s Initial Position = " << planet1.position[0] << "," <<planet1.position[1]<<","<<planet1.position[2]<< endl;
    cout << planet2.name << "'s Mass = " <<planet2.mass<< endl;
    cout << planet2.name <<"'s Initial Position = " << planet2.position[0] << "," <<planet2.position[1]<<","<<planet2.position[2]<< endl;

    //Setup the binary system
    ODEsolver binary;         //create object of class ODEsolver with default constructor ODEsolver(). If put the () in declaration here, it doesn't work!
    binary.add(planet1);      //add planets to the solver
    binary.add(planet2);
    //Tests of the setup
    cout << "Gconst = " <<binary.Gconst << endl;
    cout << "Number of Planets = " <<binary.total_planets<<endl;
    //test planet_names (vector of strings)
    for(int i=0; i<binary.total_planets;i++){
    cout << "Planet" << i+1 << "'s name is " << binary.planet_names[i] <<endl;
    }

    //Euler method
    //start clock timer
    high_resolution_clock::time_point start1 = high_resolution_clock::now();

    binary.Euler(integration_points, final_time);  //Run Euler's method ODEsolver

    //stop clock timer and output time duration
    high_resolution_clock::time_point finish1 = high_resolution_clock::now();
    duration<double> time1 = duration_cast<duration<double>>(finish1-start1);
    cout << "Euler Solver CPU time = " << time1.count() << endl;


    //NOTE: CURRENTLY THE SUN MOVES WIERDLY. VERY SMALL, OSCILLATORY MOTIONS. Can plot seperately
    //from earth to see this

    // Velocity Verlet
    //start clock timer
    high_resolution_clock::time_point start2 = high_resolution_clock::now();

    binary.VelocityVerlet(integration_points, final_time); //Run VVerlet ODEsolver

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

