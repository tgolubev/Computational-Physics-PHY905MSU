//This is the three-body problem of the Sun, Earth, and Jupiter. The Sun is still kept fixed as the center of mass of system (and at the origin of our
//coordinate system). The Velocity Verlet method is used. (Euler method can also be used by uncommenting the Euler section).
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
    int integration_points = 100000;  // # of integration points (10,000 takes a few seconds, 100,000 takes ~30sec)
    double final_time = 50.;          // End time of calculation, end time affects energy conservation
    bool corrections = false;         //Apply relativistic corrections?
    bool sun_fixed =true;

    //Set-up planets, Using April 1st, positions and velocities in Au/year. Using Sun as origin.
    //NOTE: The Sun is at origin of coordinate system.
    planet planet1("Sun",1.,0.,0.,0.,0.,0.,0.);              // planet1 (name,mass,x,y,z,vx,vy,vz), name must be in " " marks
    planet planet2("Earth",0.000003,-9.801419762323537e-01, -1.944453400829474e-01,  1.558598318609842e-05,1.12289, -6.18859, 0.000240828);
    planet planet3("Jupiter",1000*0.00095,-5.202455845386353, -1.639519984377934,  0.1232197498880794, 0.795862, -2.50117, -0.00743178);

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
           << three_body.all_planets[i].position[1]<<","<< three_body.all_planets[i].position[2]<< endl<<endl;
    }
    //Tests of the setup
    //cout << "Gconst = " <<three_body.Gconst << endl;
    //cout << "Number of Planets = " <<three_body.total_planets<<endl;
    /*test planet_names (vector of strings)
    //for(int i=0; i<binary.total_planets;i++){
    //cout << "Planet" << i+1 << "'s name is " << binary.planet_names[i] <<endl;
    }*/

    /*
    //Euler method
    //start clock timer
    high_resolution_clock::time_point start1 = high_resolution_clock::now();

    three_body.Euler(integration_points, final_time, corrections);  //Run Euler's method ODEsolver

    //stop clock timer and output time duration
    high_resolution_clock::time_point finish1 = high_resolution_clock::now();
    duration<double> time1 = duration_cast<duration<double>>(finish1-start1);
    cout << "Euler Solver CPU time = " << time1.count() << endl;
    */

    // Velocity Verlet
    //start clock timer
    high_resolution_clock::time_point start2 = high_resolution_clock::now();

    three_body.VelocityVerlet(integration_points, final_time, corrections, sun_fixed); //Run VVerlet ODEsolver
    //check if sun fixed
    cout<<"Final sun position" << three_body.all_planets[0].position[0];
    cout<<three_body.all_planets[0].position[1]<<three_body.all_planets[0].position[2]<<endl;

    //stop clock timer and output time duration
    high_resolution_clock::time_point finish2 = high_resolution_clock::now();
    duration<double> time2 = duration_cast<duration<double>>(finish2-start2);
    cout << "Velocity Verlet Solver CPU time = " << time2.count() << endl;

    return 0;
}

