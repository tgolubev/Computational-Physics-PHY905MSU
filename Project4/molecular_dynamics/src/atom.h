#ifndef ATOM_H
#define ATOM_H
#include "math/vec3.h"

class Atom
{
private:
    double m_mass;
    vec3 m_initial_position;
public:
    vec3 position;
    vec3 velocity;
    vec3 force;
    vec3 num_bndry_crossings;  //this is to keep track of boundary crossings for calculating diffusion coeff.

    Atom(double mass);
    void setInitialPosition(double x, double y, double z);
    void resetForce();
    void resetVelocityMaxwellian(double temperature);
    double mass() { return m_mass; }
    void setMass(double mass) { m_mass = mass; }
    double initial_position(int j){return m_initial_position[j];}
};
#endif
