#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <random>
#include <chrono>
#include <time.h>
//#include "planet.h"   odesolver.h already includes planet.h
#include "odesolver.h"

using namespace std;


void print_initial(int dimension, double time_step, double final_time, double *x_initial, double *v_initial, int N);
void print_final(int dimension, double *x_final, double *v_final);

int main()
{    
    // Numerical setup
    int integration_points = 10000;  // No. of integration points
    double final_time = 50.;       // End time of calculation
    //int dimension;           // No. of spatial dimensions

    //double x[3],v[3];

    //Set-up planets
    planet earth(0.000003,1.,0.0,0.0,0.0,6.3,0.); // Earth: (mass,x,y,z,vx,vy,vz)
    planet sun(1.,0.,0.,0.,0.,0.,0.);           // Sun: (mass,x,y,z,vx,vy,vz)

    //check that values from earth were actually passed to the object and can read from object. Object properties declared in planet.h
    //and defined in planet.cpp
    cout << "Earth Mass = " <<earth.mass<< endl;
    cout << "Earth Initial Position = " << earth.position[0] << "," <<earth.position[1]<<","<<earth.position[2]<< endl;

    //Euler method
    ODEsolver earth_sun(0);   //create object of class ODEsolver. The 0 is just to get constructor to work.
    earth_sun.add(sun);      //add planets to the solver
    earth_sun.add(earth);

    //tests
    cout << "Gconst = " <<earth_sun.Gconst << endl;
    cout << "Number of Planets = " <<earth_sun.total_planets<<endl;

    earth_sun.Euler(integration_points, final_time);  //Apply Euler's method

/*  // RK4
    solver binary_rk(5.0);
    binary_rk.add(planet1);
    binary_rk.add(planet2);

    for(int j=0;j<dimension;j++){
         x[j] = planet1.position[j];
         v[j] = planet1.velocity[j];
     }
 */

    /*
    // VV
    solver binary_vv(5.0);     //create object called binary_vv, member of the class solver, pass to it the radi (radius OF THE PLANET)
                               //this is  the  solver for binary (2 planet) system
    binary_vv.add(planet1);
    binary_vv.add(planet2);

    print_initial(dimension,time_step,final_time,x,v,integration_points);

    // Evolution of binary system
    /*
    //RK4
    cout << endl << "RK4: " << endl;
    binary_rk.RungeKutta4(dimension,integration_points,final_time,force,simple,1,0.);

    for(int j=0;j<dimension;j++){
        x[j] = binary_rk.all_planets[0].position[j];
        v[j] = binary_rk.all_planets[0].velocity[j];
    }
    print_final(dimension,x,v);

    cout << endl << "VV:" << endl;
    binary_vv.VelocityVerlet(dimension,integration_points,final_time,force,simple,1,0.);

    for(int j=0;j<dimension;j++){
        x[j] = binary_vv.all_planets[0].position[j];    //sets x array values equal to positions of 1st planet (all_planets[0])
        v[j] = binary_vv.all_planets[0].velocity[j];
    }
    print_final(dimension,x,v);
    */






    return 0;
}
//---------------------------------------------------------------------------------------------------------------
//Print functions

void print_initial(int dimension,double time_step, double final_time,double *x_initial,double *v_initial, int N){
    // A function that prints out the set up of the calculation

    cout << "Time step = " << time_step << "; final time = " << final_time << "; integration points = " << N << endl;

    cout << "Initial position = ";
    for(int j=0;j<dimension;j++) cout << x_initial[j] << " ";
    cout << endl;

    cout << "Initial velocity = ";
    for(int j=0;j<dimension;j++) cout << v_initial[j] << " ";
    cout << endl;
}

void print_final(int dimension,double *x_final,double *v_final){
    // A function that prints out the final results of the calculation

    cout << "Final position = ";
    for(int j=0; j<dimension; j++) cout << x_final[j] << " ";
    cout << endl;

    cout << "Final velocity = ";
    for(int j=0; j<dimension; j++) cout << v_final[j] << " ";
    cout << endl;
}


