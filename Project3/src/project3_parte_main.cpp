//This is the three-body problem of the Sun, Earth, and Jupiter. The Sun is still kept fixed as the center of mass of system (and at the origin of our
//coordinate system).


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
    int integration_points = 10000;  // No. of integration points (10,000 takes a few seconds, 100,000 takes ~30sec)
    double final_time = 50.;       // End time of calculation

    //Set-up planets
    //NOTE: right now the Sun is at origin of coordinate system, later will set the true COM of solar system to be the origin.
    planet planet1("Sun",1.,0.,0.,0.,0.,0.,0.);              // planet1 (name,mass,x,y,z,vx,vy,vz), name must be in " " marks
    planet planet2("Earth",0.000003,1.,0.,0.,0.,6.3,0.);
    planet planet3("Jupiter",0.00095,5.2,0.,0,0.,-6.3,0.);  //VELOCITY NEEDs TO BE CHANGED TO REFLECT REAL JUPITER!

    //Setup the three body system
    ODEsolver three_body;         //create object of class ODEsolver with default constructor ODEsolver(). If put the () in declaration here, it doesn't work!
    three_body.add(planet1);      //add planets to the solver
    three_body.add(planet2);
    three_body.add(planet3);

    //Output the properties of the planets
    for(int i=0;i<three_body.total_planets;i++)
    {
      cout << three_body.all_planets[i].name << "'s Mass = " <<three_body.all_planets[i].mass<< endl;
      cout << three_body.all_planets[i].name <<"'s Initial Position = " << three_body.all_planets[i].position[0] << ","
           << three_body.all_planets[i].position[1]<<","<< three_body.all_planets[i].position[2]<< endl;
    }
    //Tests of the setup
    cout << "Gconst = " <<three_body.Gconst << endl;
    cout << "Number of Planets = " <<three_body.total_planets<<endl;

    //test planet_names (vector of strings)
    //for(int i=0; i<binary.total_planets;i++){
    //cout << "Planet" << i+1 << "'s name is " << binary.planet_names[i] <<endl;
    //}

    //Euler method
    //start clock timer
    high_resolution_clock::time_point start1 = high_resolution_clock::now();

    three_body.Euler(integration_points, final_time);  //Run Euler's method ODEsolver

    //stop clock timer and output time duration
    high_resolution_clock::time_point finish1 = high_resolution_clock::now();
    duration<double> time1 = duration_cast<duration<double>>(finish1-start1);
    cout << "Euler Solver CPU time = " << time1.count() << endl;

    //NOTE: CURRENTLY THE SUN MOVES. When add jupiter, it moves in a line, quite significantly.

    // Velocity Verlet
    //start clock timer
    high_resolution_clock::time_point start2 = high_resolution_clock::now();

    three_body.VelocityVerlet(integration_points, final_time); //Run VVerlet ODEsolver

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

