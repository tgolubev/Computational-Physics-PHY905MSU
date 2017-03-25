#include "planet.h"

planet::planet()
{
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

planet::planet(double M, double x, double y, double z, double vx, double vy, double vz)
{
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

double planet::GravitationalForce(planet otherPlanet,double Gconst)
{
    double r = this->distance(otherPlanet);
    if(r!=0) return Gconst*this->mass*otherPlanet.mass/(r*r);   //F = Gm1m2/r^2
    else return 0;
}

double planet::X_GravitationalForce(planet otherPlanet,double Gconst)
{
    double r = this->distance(otherPlanet);
    double relative_x = otherPlanet.position[0]-this->position[0];
    if(r!=0) return Gconst*this->mass*otherPlanet.mass*relative_x/(r*r*r);   //F = Gm1m2/r^2
    else return 0;
}

double planet::Y_GravitationalForce(planet otherPlanet,double Gconst)
{
    double r = this->distance(otherPlanet);
    double relative_y = otherPlanet.position[1]-this->position[1];
    if(r!=0) return Gconst*this->mass*otherPlanet.mass*relative_y/(r*r*r);
    else return 0;
}

double planet::Z_GravitationalForce(planet otherPlanet,double Gconst)
{
    double r = this->distance(otherPlanet);
    double relative_z = otherPlanet.position[2]-this->position[2];
    if(r!=0) return Gconst*this->mass*otherPlanet.mass*relative_z/(r*r*r);
    else return 0;
}

double planet::Acceleration(planet otherPlanet, double Gconst)
{
    double r = this->distance(otherPlanet);
    if(r!=0) return this->GravitationalForce(otherPlanet,Gconst)/(this->mass*r);  // a = F/m1  "return this->GravitationalForce" means it will
                                                                                  //find and return gravitationalforce on the current object due to the "otherPlanet"

    else return 0;
}

double planet::KineticEnergy()
{
    double velocity_sqrd = (this->velocity[0]*this->velocity[0]) + (this->velocity[1]*this->velocity[1]) + (this->velocity[2]*this->velocity[2]);
    return 0.5*this->mass*velocity_sqrd;
}

double planet::PotentialEnergy(planet &otherPlanet, double Gconst)
{
    return -Gconst*this->mass*otherPlanet.mass/this->distance(otherPlanet);      //U = -Gm1m2/r
    //if(epsilon==0.0) return -Gconst*this->mass*otherPlanet.mass/this->distance(otherPlanet);      //U = -Gm1m2/r
    //else return (Gconst*this->mass*otherPlanet.mass/epsilon)*(atan(this->distance(otherPlanet)/epsilon) - (0.5*M_PI));  //WHAT IS THIS FORMULA??
}
