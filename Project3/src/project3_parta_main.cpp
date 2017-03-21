//Include statements for header files
#include <iostream>
#include <fstream>
#include <iomanip>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>    //  Mathematical functions like sin, cos etc
#include <string>   //  useful library to operate on characters
#include <chrono>  //for high resolution clock
//#include <armadillo.h>

using namespace std;
using namespace std::chrono; //for high res clock

ofstream output;

//void PrintInitialValues(int, double, double, double *, double *, int);
//void PrintFinalValues(int, double *, double *);

int main(int argc, char *argv[])
{
    string filename;
    int IntegrationPoints;

    if( argc <= 2 ){
         cout << "Bad Usage: "
              " Input must be: filename IntegrationPoints" << endl;
          exit(1);
     }
     else{

     filename = argv[1];
     IntegrationPoints = stoi(argv[2]); //convert to int
    }

    //int IntegrationPoints = 10000;  // No. of integration points
    double FinalTime = 50.;       // End time of calculation
    //int Dimension = 3;           // No. of spatial dimensions
    const double pi= acos(-1.0);  //define pi

    cout << "Earth-Sun binary system" << endl;

    double TimeStep = FinalTime/((double) IntegrationPoints);
    /*
    double earth_position[3],earth_velocity[3];  // positions and velocities
    earth_position[1]=1.0;
    earth_position[2]=0.0;
    earth_position[3]=0.0;

    earth_velocity[1]=0.0;
    earth_velocity[2]=2*pi;
    earth_velocity[3]=0.0;
*/
//Dynamic  allocating the arrays
    double *x = new double[IntegrationPoints+1];
    double *y = new double[IntegrationPoints+1];
    double *z = new double[IntegrationPoints+1];
    double rr;
    double rcubed;
    const double fourpi_squared = 4*pi*pi;
    double *v_x= new double[IntegrationPoints+1];
    double *v_y= new double[IntegrationPoints+1];
    double *v_z= new double[IntegrationPoints+1];

    //Initial values
    x[0]=1.0;
    y[0]=0.0;
    z[0]=0.0;
    v_x[0]=0.0;
    v_y[0]=2*pi;
    v_z[0]=0.0;


//Euler's Method
    /*
    for(int i=0; i<IntegrationPoints; i++)
    {
        x[i+1]=x[i]+TimeStep*v_x[i];
        y[i+1]=y[i]+TimeStep*v_y[i];
        z[i+1]=z[i]+TimeStep*v_z[i];
        rr = (x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
        rcubed = pow(rr,1.5);    //power function
        v_x[i+1]=v_x[i]-(fourpi_squared*x[i]*TimeStep)/rcubed;
        v_y[i+1]=v_y[i]-(fourpi_squared*y[i]*TimeStep)/rcubed;
        v_z[i+1]=v_z[i]-(fourpi_squared*z[i]*TimeStep)/rcubed;
    }\
*/


    //Velocity Verlet
    double TimeStep_sqrd = TimeStep*TimeStep;
    double ith_accel_x, ith_accel_y, ith_accel_z, next_accel_x, next_accel_y, next_accel_z;

    for(int i=0; i<IntegrationPoints; i++)
    {
        //attempt at making more efficient: use previous iteration's values except for i=0 iteration.
        /*
        if (i==0)
        {
            rr = (x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
            rcubed = pow(rr,1.5);    //power function
            ith_accel_x = -(fourpi_squared*x[i])/rcubed;
            ith_accel_y = -(fourpi_squared*y[i])/rcubed;
            ith_accel_z = -(fourpi_squared*z[i])/rcubed;
        }
        else
        {
        ith_accel_x = next_accel_x;
        ith_accel_x = next_accel_x;
        ith_accel_x = next_accel_x;
        }
        */

        rr = (x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
        rcubed = pow(rr,1.5);    //power function
        ith_accel_x = -(fourpi_squared*x[i])/rcubed;
        ith_accel_y = -(fourpi_squared*y[i])/rcubed;
        ith_accel_z = -(fourpi_squared*z[i])/rcubed;
        x[i+1]=x[i]+TimeStep*v_x[i]+0.5*TimeStep_sqrd*ith_accel_x;
        y[i+1]=y[i]+TimeStep*v_y[i]+0.5*TimeStep_sqrd*ith_accel_y;
        z[i+1]=z[i]+TimeStep*v_z[i]+0.5*TimeStep_sqrd*ith_accel_z;

        rr = (x[i+1]*x[i+1]+y[i+1]*y[i+1]+z[i+1]*z[i+1]);
        rcubed = pow(rr,1.5);    //power function
        next_accel_x = -(fourpi_squared*x[i+1])/rcubed;
        next_accel_y =-(fourpi_squared*y[i+1])/rcubed;
        next_accel_z = -(fourpi_squared*z[i+1])/rcubed;
        v_x[i+1]=v_x[i]+0.5*TimeStep*(ith_accel_x+next_accel_x);
        v_y[i+1]=v_y[i]+0.5*TimeStep*(ith_accel_y+next_accel_y);
        v_z[i+1]=v_z[i]+0.5*TimeStep*(ith_accel_z+next_accel_z);
    }


    //PrintInitialValues(Dimension,TimeStep,FinalTime,earth_position,earth_velocity,IntegrationPoints);

    cout << "Velocity Verlet results for the Sun-Earth system:" << endl;



     //PrintFinalValues(Dimension,earth_position,earth_velocity);

    //Setup output file
    output.open(filename);                                      //after initial opening of file, it's referred to by the ofstream object name.
    output << setiosflags(ios::showpoint | ios::uppercase);     //sets to write i.e. 10^6 as E6

    for (int i = 0; i <=IntegrationPoints;i++) {

       output << setw(15) << setprecision(8) << i*TimeStep;
       output << setw(15) << setprecision(8) << x[i];
       output << setw(15) << setprecision(8) << y[i];
       output << setw(15) << setprecision(8) << z[i]<<endl;
    }
    output.close();
    delete [] x; delete [] y; delete []z; delete []v_x; delete []v_y; delete []v_z; //memory that was occupied by these arrays is now freed


   return 0;
}

//use pointers to pass array into a function as in below

void PrintInitialValues(int Dimension,double TimeStep, double FinalTime,double *x_initial,double *v_initial, int N){
        // A function that prints out the set up of the calculation

        cout << "Time step = " << TimeStep << "; final time = " << FinalTime << "; integration points = " << N << endl;

        cout << "Initial position = ";
        for(int j=0;j<Dimension;j++) cout << x_initial[j] << " ";
        cout << endl;

        cout << "Initial velocity = ";
        for(int j=0;j<Dimension;j++) cout << v_initial[j] << " ";
        cout << endl;
    }

void PrintFinalValues(int Dimension,double *x_final,double *v_final){
        // A function that prints out the final results of the calculation

        cout << "Final position = ";
        for(int j=0; j<Dimension; j++) cout << x_final[j] << " ";
        cout << endl;

        cout << "Final velocity = ";
        for(int j=0; j<Dimension; j++) cout << v_final[j] << " ";
        cout << endl;
    }
