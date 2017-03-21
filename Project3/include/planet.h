//The header file decleares all the variables, constructors, destructors, and functions that are used in the class. The .cpp file for the class
//will contain the code which describes what these functions actually do.


#ifndef PLANET_H   //#ifndef directive, you can include a block of text only if a particular expression (PLANET_H in our case) is undefined;
                   //then, within the header file, you can define the expression. This ensures that the code in the header file referenced by
                   //#ifndef is included only the first time the file is loaded.


#define PLANET_H   //defines PLANET_H if it has not yet been defined (checked by #ifndef)
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
using std::vector;


class planet
{
public:
    // Properties
    double mass;
    double position[3];  //array
    double velocity[3];
    double potential;
    double kinetic;

    // Initializers (aka constructers) for the objects (planet) of the class
    planet();                                                                    //initialize by passing no values to it
    planet(double M,double x,double y,double z,double vx, double vy,double vz);  //initialize by passing these variables

    // Functions
    double distance(planet otherPlanet);
    double GravitationalForce(planet otherPlanet, double Gconst);
    double Acceleration(planet otherPlanet, double Gconst);
    double KineticEnergy();
    double PotentialEnergy(planet &otherPlanet, double Gconst, double epsilon);

};

#endif // This marks the end of the block of code which will be skipped if it was already defined (checked by #ifndef)
