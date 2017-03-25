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

int main()
{    
    // Numerical setup
    int integration_points = 10000;  // No. of integration points
    double final_time = 50.;       // End time of calculation
    //int dimension;           // No. of spatial dimensions

    //Set-up planets
    planet planet1("Earth",0.000003,1.,0.0,0.0,0.0,6.3,0.); // planet1 (name,mass,x,y,z,vx,vy,vz)   //name must be in " " marks
    planet planet2("Sun",1.,0.,0.,0.,0.,0.,0.);           // planet2 (name,mass,x,y,z,vx,vy,vz)

    //Output the properties of the planets
    cout << planet1.name << "'s Mass = " <<planet1.mass<< endl;
    cout << planet1.name <<"'s Initial Position = " << planet1.position[0] << "," <<planet1.position[1]<<","<<planet1.position[2]<< endl;
    cout << planet2.name << "'s Mass = " <<planet2.mass<< endl;
    cout << planet2.name <<"'s Initial Position = " << planet2.position[0] << "," <<planet2.position[1]<<","<<planet2.position[2]<< endl;

    //Euler method
    ODEsolver binary;         //create object of class ODEsolver with default constructor ODEsolver(). If put the () in declaration here, it doesn't work!
    binary.add(planet1);      //add planets to the solver
    binary.add(planet2);

    //tests
    cout << "Gconst = " <<binary.Gconst << endl;
    cout << "Number of Planets = " <<binary.total_planets<<endl;

    binary.Euler(integration_points, final_time);  //Apply Euler's method

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

