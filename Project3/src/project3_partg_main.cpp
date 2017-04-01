//This studies the perihelion precession of Mercury around the Sun with no other planets present. The Sun is taken to be the COM of the system.
//This makes use of the planet and ODEsolver classes.

//No input from the command line is required.

//Coded by: Tim Golubev, Hao Lin, Xingze Mao

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
    int integration_points = 100000;  // No. of integration points (10,000 takes a few seconds, 100,000 takes ~30sec). Mercury requires 500k if simulate 100 earth years for more consistent thin orbit.
    double final_time = 100.;       // End time of calculation. (Mercury takes 88days for 1 orbit)
    bool corrections = true;      //Relativistic corrections to forces

    //Set-up planets, Using April 1st, positions and velocities in Au/year. Using Sun as origin.
    //NOTE: The Sun is at origin of coordinate system.
    planet planet1("Sun",1.,0.,0.,0.,0.,0.,0.);              // planet1 (name,mass,x,y,z,vx,vy,vz), name must be in " " marks
    planet planet2("Mercury",1.66012e-7,0.3075,0.,0.,0.,12.44,0.);  //using #'s given in part g

    //Setup the system
    ODEsolver precession;         //create object of class ODEsolver with default constructor ODEsolver(). If put the () in declaration here, it doesn't work!
    precession.add(planet1);      //add planets to the solver
    precession.add(planet2);

    //Output the properties of the planets
    for(int i=0;i<precession.total_planets;i++)
    {
      cout << precession.all_planets[i].name << "'s Mass = " <<precession.all_planets[i].mass<< endl;
      cout << precession.all_planets[i].name <<"'s Initial Position = " << precession.all_planets[i].position[0] << ","
           << precession.all_planets[i].position[1]<<","<< precession.all_planets[i].position[2]<< endl;
    }
    //Tests of the setup
    cout << "Gconst = " <<precession.Gconst << endl;
    cout << "Number of Planets = " <<precession.total_planets<<endl;
    /*test planet_names (vector of strings)
    //for(int i=0; i<binary.total_planets;i++){
    //cout << "Planet" << i+1 << "'s name is " << binary.planet_names[i] <<endl;
    }*/

    /*
    //Euler method
    //start clock timer
    high_resolution_clock::time_point start1 = high_resolution_clock::now();

    precession.Euler(integration_points, final_time, corrections);  //Run Euler's method ODEsolver

    //stop clock timer and output time duration
    high_resolution_clock::time_point finish1 = high_resolution_clock::now();
    duration<double> time1 = duration_cast<duration<double>>(finish1-start1);
    cout << "Euler Solver CPU time = " << time1.count() << endl;
    */

    // Velocity Verlet
    //start clock timer
    high_resolution_clock::time_point start2 = high_resolution_clock::now();

    precession.VelocityVerlet(integration_points, final_time, corrections); //Run VVerlet ODEsolver

    //stop clock timer and output time duration
    high_resolution_clock::time_point finish2 = high_resolution_clock::now();
    duration<double> time2 = duration_cast<duration<double>>(finish2-start2);
    cout << "Velocity Verlet Solver CPU time = " << time2.count() << endl;

    return 0;
}

