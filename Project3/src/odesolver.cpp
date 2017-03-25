#include "planet.h"
#include "odesolver.h"
#include <iostream>
#include <cmath>
#include <time.h>

//NOTE: for some reason it gave error when tried to use constructor w/o passing anything
//to it. So made constructor pass zero to it

ODEsolver::ODEsolver(int zero)
{
    total_planets = zero;
    Gconst = 4*M_PI*M_PI;       //M_PI is the pi in c++: #define _USE_MATH_DEFINES
    totalKinetic = 0;
    totalPotential = 0;

}


//functions
void ODEsolver::add(planet newplanet)        //calling add, pass a planet object to it, inside the fnc add the object is named "newplanet"
{
    total_planets += 1;                   // "+=" means total_planets = total_planets +1
    all_planets.push_back(newplanet);    //push_back() is fnc. of the C++ vector class which adds a new element to the end of a vector
                                         //We wrote: #include <vector>  and using std::vector;  in planet.h so can immediately use push_back w/o ::vector
                                         //all_planets is defined as vector property of the objects of our solver class.
}


void ODEsolver::print_position(std::ofstream &output, double time,int number)
{   // Writes mass, position and velocity to a file "output"
        for(int i=0;i<number;i++){
            planet &current = all_planets[i];
            //MAYBE IF WE SET PRECISION HERE CAN GET THE OUTPUT FILES TO LOOK NICER

            output << time << "\t" << i+1 << "\t" << current.mass;         //"\t" is a tab. i.e. "Hello\tWorld" will output "Hello     World"
            for(int j=0;j<3;j++) output << "\t" << current.position[j];    //for 3D
            for(int j=0;j<3;j++) output << "\t" << current.velocity[j];
            output << std::endl;
        }
}

void ODEsolver::print_energy(std::ofstream &output, double time)
{   // Writes energies to a file "output"

    this->KineticEnergySystem();
    this->PotentialEnergySystem();
    for(int nr=0;nr<total_planets;nr++){
        planet &Current = all_planets[nr];
        output << time << "\t" << nr << "\t";
        output << Current.kinetic << "\t" << Current.potential << std::endl;
    }
}

void ODEsolver::Euler(int IntegrationPoints, double FinalTime)
{

  double TimeStep = FinalTime/((double) IntegrationPoints);
  //std::vector<planet>::const_iterator j;

  // Create files for data storage
  char *filename = new char[1000];   //set up dynamiccally alocated character string w/ pointer pointing to memory adress of "filename"
  char *filenameE = new char[1000];

  //sprintf notation: If a = 5, b = 3, then "(filename, "%d plus %d is %d", a, b, a+b);" will return "5 plus 3 is 8".
  //replaces the %d's with values of the variables that follow the string in " " in order that the variables are listed.

  //%.3f means to print the variable as a floating point decimal with a precision of at least 3 digits.
  //So i.e.: sprintf(filename, "cluster_VV_%d_%.3f.txt",total_planets,time_step); will replace the %.3f with time_step to
  //3 digits precision, so i.e. if it's 0.15, it will display as 0.150. If 0.005, it displays as 0.005

  sprintf(filename, "Euler_%d_Planets_%.3f_Positions.txt",total_planets,TimeStep);
  sprintf(filenameE, "Euler_%d_Planets_%.3f_Energies.txt",total_planets,TimeStep); //make filenameE = 'Euler_[total#ofplanets}_planets_[timestep, displayed to 3 digits].txt'

  //define objects of class ofstream which will be for the output files. Each object has a filename.
  std::ofstream output_file(filename);
  std::ofstream output_energy(filenameE);  //i.e. filenameE is what the output_energy file will be named.

  // Initialize forces and time
  double Fx,Fy,Fz; // Forces in each dimension
  double time = 0.0;

  // Write initial values to file
  print_position(output_file,time,total_planets);
  print_energy(output_energy,time);

  for(int i=0; i<IntegrationPoints; i++)   //loop for time steps
  {
    //for(j = all_planets.begin(); j!=all_planets.end(); ++j)   //for loop for vectors
    //Here to save active memory, we will not save all of the values at all time steps. We
    //will use the previous time step values to compute next time step values and output values
    //to file in each time step.

      time +=TimeStep;  //add 1 timestep each iteration, starting from i=0 iteration
      //std::cout<<time<<std::endl;

      int current_index;    //Declare outside of for loop because want to use in more than 1 loop
    for(current_index=0; current_index<total_planets; current_index++){   //loop over planets
       planet &current = all_planets[current_index];  //the & IS NECESSARY TO BE ABLE TO CHANGE THE VALUES of the object!

       Fx = Fy = Fz = 0.0; // Reset forces before each run

       for(int n=0; n<total_planets; n++){     //Calculate pairwise grav. force
           if(n==current_index)continue;                   //skip this case
           planet &other = all_planets[n];   //the & IS NECESSARY TO BE ABLE TO CHANGE THE VALUES of the object!
           Fx += current.X_GravitationalForce(other, Gconst);
           Fy += current.Y_GravitationalForce(other, Gconst);
           Fz += current.Z_GravitationalForce(other, Gconst);
       }

       std::cout << "Fx = " << Fx <<std::endl;
       //FX IS NOT CHANGING CURRENTLY!
       //update positions by 1 time step
       current.position[0]=current.position[0]+TimeStep*current.velocity[0];
       current.position[1]=current.position[1]+TimeStep*current.velocity[1];
       current.position[2]=current.position[2]+TimeStep*current.velocity[2];

       //update velocities by 1 time step
       //NOT SURE HERE IF SHOULD BE - OR +
       current.velocity[0] = current.velocity[0] - (Fx/current.mass)*TimeStep;
       current.velocity[1] = current.velocity[1] - (Fy/current.mass)*TimeStep;
       current.velocity[2] = current.velocity[2] - (Fz/current.mass)*TimeStep;

       std::cout<<"velocity = " << current.velocity[0] <<std::endl;
    }
    //print the current values to output file
    print_position(output_file,time,total_planets);
    //std::cout<<"time = " << time <<std::endl;
    print_energy(output_energy,time);

       //std::cout<<"position = " << current.position[0]<<std::endl;
    }

  // Close files
  output_file.close();
  output_energy.close();
}



void ODEsolver::KineticEnergySystem()
{
    totalKinetic = 0;
    for(int n=0;n<total_planets;n++){
        planet &Current = all_planets[n];
        Current.kinetic += Current.KineticEnergy();
    }
}

void ODEsolver::PotentialEnergySystem()
{
    totalPotential = 0;
    for(int n=0;n<total_planets;n++){         //reset all potential energies of planets to 0
        planet &Current = all_planets[n];
        Current.potential = 0;
    }
    for(int n1=0;n1<total_planets;n1++){
        planet &Current = all_planets[n1];
        for(int n2=n1+1;n2<total_planets;n2++){
            planet &Other = all_planets[n2];
            Current.potential += Current.PotentialEnergy(Other,Gconst);
            Other.potential += Other.PotentialEnergy(Current,Gconst);
        }
    }
}
