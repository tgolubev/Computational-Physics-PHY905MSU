//This is the planet class for the solar system model. It allows to create planet objects and allows to calculate the necessary
//quantities needed for modeling orbiting planets.
//The gravitational force functions include a relativistic correction which can be turned on and off by
//the bool value which is passed to the function.

//Adapted by Tim Golubev, Hao Lin, Xingze Mao from Morten Hjorth-Jensen's planet class
//found at https://github.com/CompPhysics/ComputationalPhysicsMSU


//#include <iostream>
#include "planet.h"
#include <cstring>
#include <iostream>

planet::planet()
{
    name = "nameless_planet";  //to assign hardcoded value for a string must use quotation marks
    mass = 1.;
    position[0] = 1.;
    position[1] = 0.;
    position[2] = 0.;
    velocity[0] = 0.;
    velocity[1] = 0.;
    velocity[2] = 0.;
    potential = 0.;
    kinetic = 0.;
}

planet::planet(std::string planet_name, double M, double x, double y, double z, double vx, double vy, double vz)
{
    name = planet_name;
    mass = M;
    position[0] = x;
    position[1] = y;
    position[2] = z;
    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;
    potential = 0.;
    kinetic = 0.;
}

double planet::radius()
{   //Computes the scalar radius of the planet from coordinate system origin
    double r = sqrt(this->position[0]*this->position[0]+ this->position[1]*this->position[1]+ this->position[2]*this->position[2]);
    return r;

}

double planet::distance(planet otherPlanet)   //CAREFUL!: this passes "otherPlanet" which of type/class "planet". It "planet" is not an argument!
{
    double x1,y1,z1,x2,y2,z2,xx,yy,zz;


    x1 = this->position[0];        //"this" is a constant pointer that holds the memory address of the current object ("planet")
    y1 = this->position[1];        //-> is class member access operator
    z1 = this->position[2];

    //Unless a class member name is hidden, using the class member name is equivalent to using
    //the class member name with the this pointer and the class member access operator (->).
    //Often use this-> to emphasize that accessing a member of a class.
    //another convention is to use "m_" prefix for all members of classes, i.e. m_planet and then drop the "this->"

    x2 = otherPlanet.position[0];
    y2 = otherPlanet.position[1];
    z2 = otherPlanet.position[2];

    xx = x1-x2;
    yy = y1-y2;
    zz = z1-z2;

    return sqrt(xx*xx + yy*yy + zz*zz);
 }

double planet::GravitationalForce(planet otherPlanet,double Gconst, bool relativistic)
{
    double r = this->distance(otherPlanet);
    if(r!=0 && relativistic==false) return Gconst*this->mass*otherPlanet.mass/(r*r);   //F = Gm1m2/r^2
    else if(r!=0) return (Gconst*this->mass*otherPlanet.mass*(1+this->Relativistic_correction(otherPlanet)))/(r*r);
    else return 0;
}

double planet::X_GravitationalForce(planet otherPlanet,double Gconst, bool relativistic)   //NOTE: these gravforces already have the proper sign!!
{
    double r = this->distance(otherPlanet);
    double relative_x = otherPlanet.position[0]-this->position[0];
    if(r!=0 && relativistic == false) return Gconst*this->mass*otherPlanet.mass*relative_x/(r*r*r);   //F = Gm1m2x/r^3
    else if(r!=0) return (Gconst*this->mass*otherPlanet.mass*relative_x*(1+this->Relativistic_correction(otherPlanet)))/(r*r*r);
    else return 0;
}

double planet::Y_GravitationalForce(planet otherPlanet,double Gconst, bool relativistic)
{
    double r = this->distance(otherPlanet);
    double relative_y = otherPlanet.position[1]-this->position[1];
    if(r!=0 && relativistic == false) return Gconst*this->mass*otherPlanet.mass*relative_y/(r*r*r);   //F = Gm1m2y/r^3
    else if(r!=0) return (Gconst*this->mass*otherPlanet.mass*relative_y*(1+this->Relativistic_correction(otherPlanet)))/(r*r*r);
    else return 0;
}

double planet::Z_GravitationalForce(planet otherPlanet,double Gconst, bool relativistic)
{
    double r = this->distance(otherPlanet);
    double relative_z = otherPlanet.position[2]-this->position[2];
    if(r!=0 && relativistic == false) return Gconst*this->mass*otherPlanet.mass*relative_z/(r*r*r);   //F = Gm1m2x/r^3
    else if (r!=0) return (Gconst*this->mass*otherPlanet.mass*relative_z*(1+this->Relativistic_correction(otherPlanet)))/(r*r*r);
    else return 0;
}

double planet::AngularMomentum()
{  //Calculates angular momentum using L/m = radius cross velocity, where radius is with respect to system COM
    double RcrossV_x = this->position[1]*this->velocity[2]-this->position[2]*this->velocity[1];
    double RcrossV_y = this->position[0]*this->velocity[2]-this->position[2]*this->velocity[0];
    double RcrossV_z = this->position[0]*this->velocity[1]-this->position[1]*this->velocity[0];
    double ang_mom = this->mass*sqrt(RcrossV_x*RcrossV_x + RcrossV_y*RcrossV_y + RcrossV_z*RcrossV_z);
    return ang_mom;
}

double planet::Relativistic_correction(planet otherPlanet)
{
    double const c = 63197.8;  //speed of light in: Au/year
    double AngMom_perM = this->AngularMomentum()/this->mass;
    std::cout<<"AngMom_perM = " << AngMom_perM<<std::endl;   //WHEN THESE ARE TINY, GIVES TINY  CORRECTION
    double Distance = this->distance(otherPlanet);
    std::cout<<"Distance = " <<Distance <<std::endl;
    double correction = (3*AngMom_perM*AngMom_perM)/(Distance*Distance*c*c);
    std::cout<<"correction = " <<correction << std::endl;
    return correction;
}

double planet::Acceleration(planet otherPlanet, double Gconst, bool relativistic)
{
    double r = this->distance(otherPlanet);
    if(r!=0) return this->GravitationalForce(otherPlanet,Gconst,relativistic )/(this->mass*r);
    // a = F/m1  "return this->GravitationalForce" means it will find and return gravitationalforce on the current object due to the "otherPlanet"
    else return 0;
}

double planet::Velocity_scalar()
{
    double Velocity_scalar = sqrt((this->velocity[0]*this->velocity[0]) + (this->velocity[1]*this->velocity[1]) + (this->velocity[2]*this->velocity[2]));
    return Velocity_scalar;
}

double planet::KineticEnergy()
{
    double velocity_sqrd = (this->velocity[0]*this->velocity[0]) + (this->velocity[1]*this->velocity[1]) + (this->velocity[2]*this->velocity[2]);
    return 0.5*this->mass*velocity_sqrd;
}

double planet::PotentialEnergy(planet &otherPlanet, double Gconst)
{
    return -Gconst*this->mass*otherPlanet.mass/this->distance(otherPlanet);      //U = -Gm1m2/r
}
