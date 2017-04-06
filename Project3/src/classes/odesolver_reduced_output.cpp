//This is a class of ordinary differential equation numerical solver algorithms for the solar system model.
//(The solver algorithms are universal, just the naming of variables is specific to the solar system problem.)

//Coded by: Tim Golubev, Hao Lin, Xingze Mao

#include "planet.h"
#include "odesolver.h"
#include <iostream>
#include <cmath>
#include <time.h>
#include <array>
#include <iomanip>
#include <string>       //need this to be able to cout strings


ODEsolver::ODEsolver()          //to use this constructor in code, just use i.e. : ODEsolver object_name;     //NOT: ODEsolver object_name()
{
    total_planets = 0;
    Gconst = 4*M_PI*M_PI;       //M_PI is the pi in c++: #define _USE_MATH_DEFINES
    totalKinetic = 0;
    totalPotential = 0;
}


//functions
void ODEsolver::add(planet newplanet)        //calling add, pass a planet object to it, inside the fnc add the object is named "newplanet"
{
    total_planets += 1;                      // "+=" means total_planets = total_planets +1
    all_planets.push_back(newplanet);
    planet_names.push_back(newplanet.name);
    totalKinetic += newplanet.KineticEnergy();
    /*push_back() is fnc. of the C++ vector class which adds a new element to the end of a vector
    We wrote: #include <vector>  and using std::vector;  in planet.h so can immediately use push_back w/o ::vector
    all_planets is defined as vector property of the objects of our solver class.*/
}

double **ODEsolver::save_initial_values()        //** here means it returns a pointer
{
    //Save initial values into a (total_planets, 6) size array for reset later. Columns 1-3 are initial positions, 4-6 are initial velocities
    //std::cout<<"Initial Values save test: " << std::endl;
    double **initial_values;
    initial_values = new double*[total_planets];          //each row is each planet
    for(int i=0;i<total_planets;i++){
       planet current = all_planets[i];                   //select each planet out of vector of planet objects
       initial_values[i] = new double[6];                 //3 for positions, 3 for velocities
       for(int j=0;j<6;j++){
           if(j<3){
            initial_values[i][j] = current.position[j];
            //std::cout<<initial_values[i][j]<<std::endl;
           }else{
            initial_values[i][j] = current.velocity[j-3];  //velocities are stored in elements 0-2 of velocity[]
            //std::cout<<initial_values[i][j]<<std::endl;
           }
       }
    }
    return initial_values;  //returns a pointer to initial_positions
}


void ODEsolver::reset_initial_values(double** initial_values)   //double** means are passing a pointer
{   //Reset initial values to what they were before running a solver algo
    //std::cout<<"Reset Initial values test: " <<std::endl;
    for(int i=0;i<total_planets;i++){
       planet &current = all_planets[i];      //select each planet out of vector of planet objects
       for(int j=0;j<6;j++){
           if(j<3){
            current.position[j]=initial_values[i][j];
            //std::cout<< current.position[j]<<std::endl;
           }else{
            current.velocity[j-3]=initial_values[i][j];  //velocities are stored in elements 0-2 of velocity[]
            //std::cout<< current.velocity[j-3]<<std::endl;
           }
       }
    }
}


void ODEsolver::print_position(std::ofstream &output, double time,int number)
{   // Writes mass, position and velocity to a file "output"
        for(int i=1;i<number;i++){       //Only print out the 1st planet (ie Mercury). Don't need all values for sun, which stays fixed.
            planet &current = all_planets[i];
            output << std::setiosflags(std::ios::showpoint | std::ios::uppercase);     //sets to write i.e. 10^6 as E6
            output << std::setprecision(11) << time << "\t" << i+1; //<< "\t" << current.mass;        Skip outputing mass also //"\t" is a tab. i.e. "Hello\tWorld" will output "Hello     World"
            for(int j=0;j<2;j++) output << "\t" << std::setprecision(10) << current.position[j];    //just print 2d positions, b/c z=0 for precession study
            //for(int j=0;j<3;j++) output << "\t" << std::setprecision(8) << current.velocity[j]; %To save memory, Don't print velocities since not neccesary
            output << std::endl;
        }
}

void ODEsolver::print_energy(std::ofstream &output, double time)
{   // Writes energies to a file "output"

    this->KineticEnergySystem();        //run fnc which computes KE
    this->PotentialEnergySystem();
    output << time << "\t"<<  totalKinetic << "\t" << totalPotential << "\t"<< totalKinetic+totalPotential<<std::endl;   //outputs total system energies
}

double **ODEsolver::setup_matrix(int height,int width)  //returns a double_pointer
{   // Function to set up a 2D array

    // Set up matrix
    double **matrix;
    matrix = new double*[height];

    // Allocate memory
    for(int i=0;i<height;i++)
        matrix[i] = new double[width];

    // Set values to zero
    for(int i = 0; i < height; i++){
        for(int j = 0; j < width; j++){
            matrix[i][j] = 0.0;
        }
    }
    return matrix;
}

void ODEsolver::delete_matrix(double **matrix)         //accepts a double_pointer as input
{   // Function to deallocate memory of a 2D array which  has # of rows = total_planets

    for (int i=0; i<total_planets; i++)
        delete [] matrix[i];
    delete [] matrix;
}
//--------------------------------------------------------------------------------------------------------------------------------------------
void ODEsolver::Euler(int IntegrationPoints, double FinalTime, bool corrections, bool sun_fixed)
{
  double TimeStep = FinalTime/((double) IntegrationPoints);
  double **initial = save_initial_values();              //initial contains the ADDRESS of the array of initial values
                                                         //**initial = the value of the first element of the array
  // Create files for data storage
  char *filename = new char[1000];   //set up dynamiccally alocated character string w/ pointer pointing to memory adress of "filename"
  char *filenameE = new char[1000];
  sprintf(filename, "Euler_%d_Planets_%.3f_Positions.txt",total_planets,TimeStep);
  sprintf(filenameE, "Euler_%d_Planets_%.3f_Energies.txt",total_planets,TimeStep);
  /*sprintf notation: If a = 5, b = 3, then "(filename, "%d plus %d is %d", a, b, a+b);" will return "5 plus 3 is 8".
  replaces the %d's with values of the variables that follow the string in " " in order that the variables are listed.
  %.3f means to print the variable as a floating point decimal with a precision of at least 3 digits.
  So i.e.: sprintf(filename, "Euler_%d_%.3f.txt",total_planets,time_step); will replace the %.3f with time_step to
  3 digits precision, so i.e. if it's 0.15, it will display as 0.150. If 0.005, it displays as 0.005
  */

  //define objects of class ofstream which will be for the output files. Each object has a filename.
  std::ofstream output_file(filename);
  std::ofstream output_energy(filenameE);  //i.e. filenameE is the name of the string which contains the filename

  // Initialize time and write initial values to file
  double time = 0.0;
  print_position(output_file,time,total_planets);
  print_energy(output_energy,time);

  // Initialize force components
  double Accel_x,Accel_y,Accel_z;

  //Save the initial total system energies and angular momentums
  double Initial_kinetic, Initial_potential;
  double *initial_ang_mom = new double[total_planets];
  Initial_kinetic = totalKinetic;
  std::cout<< "Initial_kinetic = " << Initial_kinetic <<std::endl;
  Initial_potential = PotentialEnergySystem();
  std::cout<< "Initial_potential = " << Initial_potential <<std::endl;
  for(int n=0;n<total_planets;n++){
        planet current = all_planets[n];
        //initial_ang_mom[n] = current.mass*current.Velocity_scalar()*current.radius();
        initial_ang_mom[n] = current.AngularMomentum();
        std::cout << all_planets[n].name << "'s Initial angular_momentum = " <<initial_ang_mom[n] << std::endl;
  }

  int current_start = 0;
  if(sun_fixed) current_start = 1;    //if sun fixed, don't time evolve the sun's position

  for(int i=0; i<IntegrationPoints; i++)   //loop over time steps
  {
    time +=TimeStep;         //add 1 timestep each iteration, starting from i=0 iteration

    int current_index;    //Declare outside of for loop over indices because want to use in more than 1 loop
    for(current_index=current_start; current_index<total_planets; current_index++){
       planet &current = all_planets[current_index];  //the & IS NECESSARY TO BE ABLE TO CHANGE THE VALUES of the object!
       Accel_x = Accel_y = Accel_z = 0.0; // Reset accelerations before each run

       for(int n=0; n<total_planets; n++){                 //Calculate pairwise grav. force
           if(n==current_index)continue;                   //skip this case
           planet &other = all_planets[n];
           Accel_x += current.X_Acceleration(other, Gconst,corrections);  //correction is a bool value (true or false).
           Accel_y += current.Y_Acceleration(other, Gconst,corrections);
           Accel_z += current.Z_Acceleration(other, Gconst,corrections);
       }
       //update positions by 1 time step
       current.position[0]+=TimeStep*current.velocity[0];
       current.position[1]+=TimeStep*current.velocity[1];
       current.position[2]+=TimeStep*current.velocity[2];

       //update velocities by 1 time step (NOTE: the forces already have the proper sign)
       current.velocity[0]+=Accel_x*TimeStep;
       current.velocity[1]+=Accel_y*TimeStep;
       current.velocity[2]+=Accel_z*TimeStep;
    }
    //print the current values to output file
    print_position(output_file,time,total_planets);
    print_energy(output_energy,time);
 }

  // Close files, reset initial values, clear memory
  output_file.close();
  output_energy.close();
  reset_initial_values(initial);
  delete_matrix(initial);  //clear memory

  //Angular Momentum conservation check
  double *ang_mom_change = new double[total_planets];
  for(int n=0;n<total_planets;n++){       //planet n=0 is center of mass of system
        planet current = all_planets[n];
        ang_mom_change[n] = current.AngularMomentum()-initial_ang_mom[n];
        if(initial_ang_mom[n]!=0.0){
        std::cout << current.name << "'s Relative angular momentum change = " <<100.0*ang_mom_change[n]/initial_ang_mom[n] << " %" << std::endl;
        }else std::cout << current.name <<"'s Initial angular momentum was 0, Angular momentum change = " << ang_mom_change[n] << std::endl;
  }
  std::cout<<std::endl;

  //Energy conservation checks
  //std::cout<< "Initial Energies = " <<Initial_kinetic << "\t" << Initial_potential <<std::endl;
  double KE_change = totalKinetic - Initial_kinetic;
  double PE_change = totalPotential -Initial_potential;
  std::cout << "Relative kinetic energy change = " << 100.0*KE_change/Initial_kinetic <<" %" << std::endl;
  std::cout <<"Relative potential energy change (Recall PE is <0) = " << 100.0*PE_change/abs(Initial_potential) << " %" <<std::endl;    //use abs b/c PE<0
  std::cout <<"Total energy loss = " << 100.0*(KE_change + PE_change)/(Initial_kinetic+Initial_potential) << " %" << std::endl<<std::endl;
  }

//---------------------------------------------------------------------------------------------------------------------------------------------------------
void ODEsolver::VelocityVerlet(int IntegrationPoints, double FinalTime, bool corrections, bool sun_fixed)
{   /*  Velocity-Verlet solver for two coupled ODEs in 3D.
    The algorithm is, exemplified in 1D for position x(t), velocity v(t) and acceleration a(t):
    x(t+dt) = x(t) + v(t)*dt + 0.5*dt*dt*a(t);
    v(t+dt) = v(t) + 0.5*dt*[a(t) + a(t+dt)];*/

    double **initial = save_initial_values();
    double TimeStep = FinalTime/((double) IntegrationPoints);
    double TimeStep_sqrd = TimeStep*TimeStep;
    double time = 0.0;

    // Create files for data storage
    char *filename = new char[1000];   //set up dynamiccally alocated character string w/ pointer pointing to memory adress of "filename"
    char *filenameE = new char[1000];
    sprintf(filename, "VVerlet_%d_Planets_%.3f_Positions.txt",total_planets,TimeStep);
    sprintf(filenameE, "VVerlet_%d_Planets_%.3f_Energies.txt",total_planets,TimeStep);
    //define objects of class ofstream which will be for the output files. Each object has a filename.
    std::ofstream output_file(filename);
    std::ofstream output_energy(filenameE);

    // Initialize forces and acceleration
    //double ith_Fx,ith_Fy,ith_Fz, next_Fx, next_Fy, next_Fz;
    double **ith_accel = setup_matrix(total_planets,3);
    //double ith_accel[3],
    double next_accel[3];

    // Write initial values to file
    //print_position(output_file,time, total_planets);
    //print_energy(output_energy,time);

    //Save the initial total system energies and angular momentums
    double Initial_kinetic, Initial_potential;
    double *initial_ang_mom = new double[total_planets];
    Initial_kinetic = totalKinetic;   //was calculated when add planets
    std::cout<< "Initial_kinetic = " << Initial_kinetic <<std::endl;
    Initial_potential = PotentialEnergySystem();
    std::cout<< "Initial_potential = " << Initial_potential <<std::endl;
    for(int n=0;n<total_planets;n++){
          planet current = all_planets[n];
          initial_ang_mom[n] = current.AngularMomentum();
          //initial_ang_mom[n] = current.mass*current.Velocity_scalar()*current.radius();
          std::cout << all_planets[n].name << "'s Initial angular_momentum = " <<initial_ang_mom[n] << std::endl;
    }

    int current_start = 0;
    if(sun_fixed) current_start = 1;    //if sun fixed, don't time evolve the sun's position

    for(int i=0; i<IntegrationPoints; i++)   //loop for time steps
    { //will use the previous time step values to compute next time step values and output values
      //to file in each time step.
      time +=TimeStep;
      int current_index;    //Declare outside of for loop over current_index because want to use in more than 1 loop

      for(current_index=current_start; current_index<total_planets; current_index++){
         planet &current = all_planets[current_index];                      //the & IS NECESSARY TO BE ABLE TO CHANGE THE VALUES of the object!
         //ith_Fx = ith_Fy = ith_Fz = next_Fx = next_Fy = next_Fz = 0.0;      // Reset forces for each time iteration and planet
         ith_accel[current_index][0]=ith_accel[current_index][1]=ith_accel[current_index][2] = 0.0;  //Reset accelerations for each time iteration and planet

         for(int n=0; n<total_planets; n++){     //Calculate pairwise grav. force
             if(n==current_index)continue;       //skip this case
             planet &other = all_planets[n];
             ith_accel[current_index][0] += current.X_Acceleration(other, Gconst, corrections);
             ith_accel[current_index][1] += current.Y_Acceleration(other, Gconst, corrections);
             ith_accel[current_index][2] += current.Z_Acceleration(other, Gconst, corrections);
         }
         //Note: the accelerations already have the proper sign!

         // Calculate new position components for current planet
         for(int j=0; j<3; j++) {
             current.position[j] += current.velocity[j]*TimeStep + 0.5*TimeStep_sqrd*ith_accel[current_index][j];
         }
      }

      //Now that have found new positions for all planets, recalculate the forces
      for(current_index=current_start; current_index<total_planets; current_index++){
         planet &current = all_planets[current_index];
         //next_Fx = next_Fy = next_Fz = 0.0;      // Reset forces for each time iteration and planet
         next_accel[0] = next_accel[1] = next_accel[2] = 0.0;   //Reset accelerations for each time iteration and planet

         for(int n=0; n<total_planets; n++){     //Calculate pairwise grav. force
             if(n==current_index)continue;       //skip this case
             planet &other = all_planets[n];
             next_accel[0] += current.X_Acceleration(other, Gconst, corrections);
             next_accel[1] += current.Y_Acceleration(other, Gconst, corrections);
             next_accel[2] += current.Z_Acceleration(other, Gconst, corrections);
         }

         // Calculate new velocity componentsfor current planet
         for(int j=0; j<3; j++) current.velocity[j] += 0.5*TimeStep*(ith_accel[current_index][j] + next_accel[j]);
     }

      if(i>0.997500000*IntegrationPoints){          //output last 0.5% of positions
      //print the current values to output file
      print_position(output_file,time,total_planets);
      //print_energy(output_energy,time);
      }
  }

  //Don't do these checks to save CPU
  /*
  //Angular Momentum conservation check
  double *ang_mom_change = new double[total_planets];
  for(int n=0;n<total_planets;n++){
        planet current = all_planets[n];
        ang_mom_change[n] = current.AngularMomentum()-initial_ang_mom[n];
        if(initial_ang_mom[n]!=0.0){
        std::cout << current.name << "'s Relative angular momentum change = " <<100.0*ang_mom_change[n]/initial_ang_mom[n] << " %" << std::endl;
        }else std::cout << current.name <<"'s Initial angular momentum was 0, Angular momentum change = " << ang_mom_change[n] << std::endl;
  }//SHOULD WE CALCULATE ANG MOMENTUM FOR ALL PLANETS TOTAL?
  std::cout<<std::endl;

  //Energy conservation checks
  //std::cout<< "Initial Energies = " <<Initial_kinetic << "\t" << Initial_potential <<std::endl;
  double KE_change = totalKinetic - Initial_kinetic;
  double PE_change = totalPotential -Initial_potential;

  std::cout << "Relative kinetic energy change = " << 100.0*KE_change/Initial_kinetic << " %" <<std::endl;
  std::cout <<"Relative potential energy change (Recall PE is <0) = " << 100.0*PE_change/abs(Initial_potential) << " %" <<std::endl;
  std::cout <<"Total energy change = " << 100.0*(KE_change + PE_change)/(Initial_kinetic+Initial_potential) << " %" << std::endl <<std::endl;
  */

  // Close files, reset initial values, clear memory
  output_file.close();
  output_energy.close();
  reset_initial_values(initial);
  delete_matrix(ith_accel);
  delete_matrix(initial);
}

double ODEsolver::KineticEnergySystem()
{
    totalKinetic = 0;
    for(int n=0;n<total_planets;n++){
        planet &Current = all_planets[n];
        Current.kinetic = Current.KineticEnergy();   //this updates the kinetic energy property of the curent planet object = 0.5mv^2
                                                     //KineticEnergy() is fnc. in planet class which calculates KE
        totalKinetic += Current.kinetic;             //adds KE of current planet to total system KE
    }
    return totalKinetic;
}

double ODEsolver::PotentialEnergySystem()
{
    totalPotential = 0;
    for(int n=0;n<total_planets;n++){         //reset all potential energies of planets to 0
        planet &Current = all_planets[n];
        Current.potential = 0;
    }
    for(int n1=0;n1<total_planets;n1++){
        planet &Current = all_planets[n1];
        for(int n2=n1+1;n2<total_planets;n2++){   //to avoid double counting
            planet &Other = all_planets[n2];
            Current.potential += Current.PotentialEnergy(Other,Gconst);
        }
        totalPotential +=Current.potential;   //add PE of current planet to total system PE.
    }
    return totalPotential;
}
