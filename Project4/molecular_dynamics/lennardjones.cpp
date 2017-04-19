#include "lennardjones.h"
#include "system.h"

double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJones::calculateForces(System &system)  //object system is passed by reference (allows changing)
{
 int current_index = 0;
 for(Atom *current_atom : system.atoms()) {
     //loop over the entire vector m_atoms (atoms() returns m_atoms: vector of pointers to the atoms objects

    for(Atom *other_atom : system.atoms()){     //Calculate pairwise forces
        if(other_atom->position==current_atom->position)continue;     //skip case of same atom
        //current_atom->force +=

    }
    current_index++;  //counts which atom are currently considering
 }
    m_potentialEnergy = 0; // Remember to compute this in the loop
}
