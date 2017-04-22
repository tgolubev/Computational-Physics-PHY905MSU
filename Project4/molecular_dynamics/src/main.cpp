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
    vec3 numberOfUnitCellsEachDimension(3,3,3);
    double initialTemperature = UnitConverter::temperatureFromSI(200.0); // measured in Kelvin
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms

    // If a first argument is provided, it is the number of unit cells and use same # of unit cells for all dimensions
    if(numberOfArguments > 1) numberOfUnitCellsEachDimension[0] = numberOfUnitCellsEachDimension[1] = numberOfUnitCellsEachDimension[2]= atoi(argumentList[1]);
    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    if(numberOfArguments > 3) latticeConstant = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));

    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds (1fs is common)

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;

    System system;
    system.createFCCLattice(numberOfUnitCellsEachDimension, latticeConstant, initialTemperature);
    //set the potential parameters
    system.potential().setEpsilon(1.0);    //i.e. LJ depth
    system.potential().setSigma(3.405);      //i.e. LJ atom diameter. Ar = 3.405 Angstroms


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

    for(int timestep=0; timestep<50000; timestep++) {  //chose # of timesteps here
        system.step(dt);   //advance system by 1 step. NOTE: PBCs ARE APPLIED IN THIS STEP: CALLS INTEGRATE WHICH IS IN velocityverlet.cpp
        //statisticsSampler.sample(system);   //use sampler to calculate system parameters
        if(timestep % 10 ==0){
            //to save CPU, don't sample every timestep
            statisticsSampler.sample(system);   //use sampler to calculate system parameters
        }
            /*
            //periodically couple system with heat bath
            //this causes wierd effect of freezing the gas!! when have a gas
            //PERHAPS INSTEAD OF RESETTING THE DISTRIBUTION EACH TIME: JUST RESCALE THE VELOCITIES TO GET THE DESIRED
            //TEMPERATURE PERIODICALLY
            for(Atom *atom : system.atoms()) {
                atom->resetVelocityMaxwellian(initialTemperature);
            }
            system.removeTotalMomentum();
            */


         if( timestep % 1000 == 0 ) {
            // Print the timestep and system properties every 1000 timesteps
            cout << setw(20) << system.steps() <<
                    setw(20) << system.time() <<
                    setw(20) << statisticsSampler.temperature() <<
                    setw(20) << statisticsSampler.kineticEnergy() <<
                    setw(20) << statisticsSampler.potentialEnergy() <<
                    setw(20) << statisticsSampler.totalEnergy() << endl;
        }
        movie.saveState(system);  //calls saveState fnc in io.cpp which saves the state to the movie.xyz file
    }

    //stop clock timer and output time duration
    high_resolution_clock::time_point finish2 = high_resolution_clock::now();
    duration<double> time2 = duration_cast<duration<double>>(finish2-start2);
    cout << "CPU time = " << time2.count() << endl;

    movie.close();

    return 0;
}
