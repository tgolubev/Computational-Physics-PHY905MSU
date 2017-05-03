//This is an object oriented molecular dynamics code which relies on many classes. This file allows to control most of the
//program parameters and choose which simulation approach (NVT, NVE, or gradual heating) to use. Certain blocks need
//to be commented/uncommented to use different approaches. This file can take input arguments, but they are not required.

//Outputs are a movie.xyz file with atom positions which can allow to create animations in i.e Ovito and statistics.txt
//which calculates system energetics, temperature, and diffusion coefficient. The frequency of output can be varied as
//desired. Of course, less frequent output, makes the program run faster.

//Author: Timofey Golubev based on the code found at https://github.com/andeplane/molecular-dynamics-fys3150.


#include "math/random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <time.h>

using namespace std;
using namespace chrono;

int main(int numberOfArguments, char **argumentList)
{
    vec3 numberOfUnitCellsEachDimension(5,5,5);
    double initialTemperature = UnitConverter::temperatureFromSI(610.0); // measured in Kelvin
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.256); // measured in angstroms

    // If a first argument is provided, it is the number of unit cells and use same # of unit cells for all dimensions
    if(numberOfArguments > 1) numberOfUnitCellsEachDimension[0] = numberOfUnitCellsEachDimension[1] = numberOfUnitCellsEachDimension[2]= atoi(argumentList[1]);
    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    if(numberOfArguments > 3) latticeConstant = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));

    double dt = UnitConverter::timeFromSI(2e-14); // Measured in seconds (1fs is common). ANYTHING LARGER THAN 2E-14 HAS ISSUES: T at step 1 is already overshooted

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;

    System system;
    system.createFCCLattice(numberOfUnitCellsEachDimension, latticeConstant, initialTemperature);
    system.potential().setEpsilon(1.0);
    system.potential().setSigma(UnitConverter::lengthFromAngstroms(3.405));      //i.e. LJ atom diameter, Ar = 3.405 Angstroms
    system.m_sample_freq=100; //statistics sampler freq.
    system.removeTotalMomentum();

    StatisticsSampler statisticsSampler;
    IO movie("movie.xyz"); // To write the state to file, create IO object called "movie"

    cout << setw(20) << "Timestep" <<
            setw(20) << "Time" <<
            setw(20) << "Temperature(not K!)" <<
            setw(20) << "KineticEnergy" <<
            setw(20) << "PotentialEnergy" <<
            setw(20) << "TotalEnergy" << endl;

    high_resolution_clock::time_point start2 = high_resolution_clock::now();  //start clock timer

    for(int timestep=0; timestep<50000000; timestep++) {  //chose # of timesteps here
        system.step(dt);   //advance system by 1 step. NOTE: PBCs ARE APPLIED IN THIS STEP: CALLS INTEGRATE WHICH IS IN velocityverlet.cpp

        //Uncomment the below block to implement gradual heating of the system.
        /*
        //heat system gradually
        statisticsSampler.sampleKineticEnergy(system);      //can't sample temperature w/o sampling KE!
        statisticsSampler.sampleTemperature(system);  //sample temp at every timestep
        system.increaseTemperature(statisticsSampler, UnitConverter::temperatureFromSI(0.0001));  //Increase T by this increment (in K)
        */

        //use sampler to calculate system parameters
        if(timestep % system.m_sample_freq ==0){
            //to save CPU, don't sample every timestep
            statisticsSampler.sample(system);
        }

        //Uncoment the below block to use NVT ensemble.
        /*
        //periodically rescale Velocities to keep T constant (NVT ensemble)
        if(timestep % 100 == 0){
            //CAN'T RESCALE MORE FREQUENTLY THAN STAT SAMPLING RATE!
           system.rescaleVelocities(statisticsSampler, initialTemperature);
        }
        */

         if( timestep % 10000 == 0 ) {
            // Print the timestep and system properties every 1000 timesteps
            cout << setw(20) << system.steps() <<
                    setw(20) << system.time() <<
                    setw(20) << statisticsSampler.temperature() <<
                    setw(20) << statisticsSampler.kineticEnergy() <<
                    setw(20) << statisticsSampler.potentialEnergy() <<
                    setw(20) << statisticsSampler.totalEnergy() << endl;
        }
        if(timestep % 1000 ==0){
          //save atom coordinates only periodically to save CPU and file size
          movie.saveState(system);  //calls saveState fnc in io.cpp which saves the state to the movie.xyz file
        }
    }

    //stop clock timer and output time duration
    high_resolution_clock::time_point finish2 = high_resolution_clock::now();
    duration<double> time2 = duration_cast<duration<double>>(finish2-start2);
    cout << "CPU time = " << time2.count() << endl;

    movie.close();

    return 0;
}
