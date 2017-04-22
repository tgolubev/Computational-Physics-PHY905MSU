#ifndef ATOM_H
#define ATOM_H
#include "math/vec3.h"

class Atom
{
private:
    double m_mass;  //NOTE: THIS WAS A FLOAT IN ORIGINAL VERSION
    vec3 initial_position;
public:
    vec3 position;
    vec3 velocity;
    vec3 force;

    Atom(double mass);
    void Atom::setInitialPosition(double x, double y, double z);
    void resetForce();
    void resetVelocityMaxwellian(double temperature);
    double mass() { return m_mass; }
    void setMass(double mass) { m_mass = mass; }
};
#endif
