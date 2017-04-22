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

double LennardJones::twntyfour_epsilon() const
{
    return 24.0*m_epsilon;
}


void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJones::calculateForces(System &system)  //object system is passed by reference (allows changing)
{
 for(Atom *current_atom : system.atoms()) {
     //loop over the entire vector m_atoms (atoms() returns m_atoms: vector of pointers to the atoms objects

     //Note: the forces for all atoms are already reset in system.cpp calculateForces() fnc.
     //check if forces are already reset
     //std::cout << "current_atom force = " <<  current_atom->force[0] << current_atom->force[1]<< current_atom->force[2]<<std::endl;
     //std::cout << current_atom <<std::endl;    //shows that current_atom is a pointer to the location of the atom position

    for(Atom *other_atom : system.atoms()){     //Calculate pairwise forces
        if(other_atom == current_atom) continue; //skip case of same atom

        //distance and vector btw objects dx should obey min. image convention: only closest distance to particle or its image is considered
        vec3 displacement;
        for(int j=0;j<3;j++){
             displacement[j] = current_atom->position[j] - other_atom->position[j];  //WHEN USE POINTERS FOR ATOMS MUST USE -> TO ACCESS THEIR PROP.'S!
             //for cases where the folded back particle will be closer than its image to a given particle
             if (displacement[j] >  system.systemSize(j) * 0.5) displacement[j] -= system.systemSize(j);   //systemSize(j) returns m_systemSize[j] from system class
             if (displacement[j] <= -system.systemSize(j) * 0.5) displacement[j] += system.systemSize(j);
         }

        //precalculate
        //declare the variables inside loop since they are only needed within the loop
        double radius = displacement.length();
        double sigma_over_radius = m_sigma/radius;
        double sigma_over_radius_sqrd = m_sigma/(radius*radius);
        double total_force = twntyfour_epsilon()*(2.0*pow(sigma_over_radius,11.)-pow(sigma_over_radius,5.))*(sigma_over_radius_sqrd); //for a single pair

        //ATTRACTIVE FORCE SHOULD POINT TOWARDS OTHER ATOM. REPULSIVE AWAY FROM OTHER ATOM!!!

        //CAN PRECALCULATE THE SIGMA POWERS: WILL PROBABLY BE MORE EFFICIENT!

        //find and set force components
        for(int j=0;j<3;j++) current_atom->force[j] += total_force*displacement[j]/radius; //i.e. Fx = F*x/r
    }

 }
    m_potentialEnergy = 0; // Remember to compute this in the loop
}
