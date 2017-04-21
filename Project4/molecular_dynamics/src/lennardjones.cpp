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
    return 24*m_epsilon;
}


void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJones::calculateForces(System &system)  //object system is passed by reference (allows changing)
{
 for(Atom *current_atom : system.atoms()) {
     //loop over the entire vector m_atoms (atoms() returns m_atoms: vector of pointers to the atoms objects
     //reset the forces
     //current_atom->force.set(0,0,0);  //using vec3 fnc set
     current_atom->force[0]=current_atom->force[1]=current_atom->force[2] = 0.0;

     //Note: the forces for all atoms are already reset in system.cpp calculateForces() fnc.

    for(Atom *other_atom : system.atoms()){     //Calculate pairwise forces
        if(other_atom == current_atom) continue;
        //if(other_atom->position==current_atom->position)continue;     //skip case of same atom

        //Note: x(), y(),z() are fncs in vec3 class which return components.
        //declare the variables inside loop since they are only needed within the loop
        //double relative_x = other_atom->position.x()-current_atom->position.x();
        //double relative_y = other_atom->position.y()-current_atom->position.y();
        //double relative_z = other_atom->position.z()-current_atom->position.z();
        //double radius = sqrt(relative_x*relative_x + relative_y*relative_y + relative_z*relative_z);

        //distance and vector btw objects dx should obey min. image convention: only closest distance to particle or its image is considered

        vec3 displacement;
        for(int j=0;j<3;j++){
             displacement[j] = other_atom->position[j] - current_atom->position[j];  //WHEN USE POINTERS FOR ATOMS MUST USE -> TO ACCESS THEIR PROP.'S!
             //for cases where the folded back particle will be closer than its image to a given particle
             if (displacement[j] >  system.systemSize(j) * 0.5) displacement[j] -= system.systemSize(j);   //systemSize(j) returns m_systemSize[j] from system class
             if (displacement[j] < -system.systemSize(j) * 0.5) displacement[j] += system.systemSize(j);
         }


        //precalculate
        double radius = displacement.length();

        double sigma_over_radius = m_sigma/radius;
        double sigma_over_radius_sqrd = m_sigma/(radius*radius);
        double force_magnitude = twntyfour_epsilon()*(2*pow(sigma_over_radius,11.)-pow(sigma_over_radius,5.))*(sigma_over_radius_sqrd); //for a single pair

        //find and set force components
        for(int j=0;j<3;j++) current_atom->force[j] += force_magnitude*displacement[j]/radius; //i.e. Fx = F*x/r
        //current_atom->force[0] += force_magnitude*relative_x/radius;
        //current_atom->force[1] += force_magnitude*relative_y/radius;
        //current_atom->force[2] += force_magnitude*relative_z/radius;

    }//THE CURRENT ISSUE IS NOT HERE: I RETESTED THIS PART!!

 }
    m_potentialEnergy = 0; // Remember to compute this in the loop
}
