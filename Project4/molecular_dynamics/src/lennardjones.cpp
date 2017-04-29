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
    m_four_epsilon = 4.0*m_epsilon;
    m_twntyfour_epsilon = 24.0*m_epsilon;  //must reset this too
}

void LennardJones::calculateForces(System &system)  //object system is passed by reference (allows changing)
{
 m_potentialEnergy = 0;  //reset potential energy

 for(int current_index=0; current_index<system.num_atoms()-1; current_index++){  //-1 b/c don't need to calculate pairs when get to last atom
    Atom *current_atom = system.atoms(current_index);  //system.atoms(index) returns the pointer to the atom corresponding to index

    for(int other_index=current_index+1;other_index<system.num_atoms();other_index++){   //to avoid double counting
        Atom *other_atom = system.atoms(other_index);
        //if(other_atom == current_atom) continue; //not necessary here, b/c it's already avoided with the for loops

        //distance and vector btw objects dx should obey min. image convention: only closest distance to particle or its image is considered
        vec3 displacement(0.,0.,0.);


        for(int j=0;j<3;j++){
             displacement[j] = current_atom->position[j] - other_atom->position[j];
             //for cases where the folded back particle will be closer than its image to a given particle
             if (displacement[j] >  system.halfsystemSize(j)) displacement[j] -= system.systemSize(j);   //systemSize(j) returns m_systemSize[j] from system class
             if (displacement[j] <= -system.halfsystemSize(j)) displacement[j] += system.systemSize(j);
         }


        //precalculate
        //declare the variables inside loop since they are only needed within the loop
        double radiusSqrd = displacement.lengthSquared();
        double radius = sqrt(radiusSqrd);

        //apply interaction range cutoff
        //std::cout<<radius<<std::endl;
        //if(radius > 5.0*m_sigma) continue;
        //THIS DOESN'T WORK FOR SOME REASON!


        //double sigma_over_radius = m_sigma/radius;
        //double sigma_over_radius_sqrd = m_sigma/radiusSqrd;
        //double sigma_over_radius_sqrd = m_sigma/(radius*radius);  //IF REPLACE ABOVE LINE WITH THIS ONE, GET DIFFERENT RESULTS!!!!!
        //double total_force = m_twntyfour_epsilon*(2.0*pow(sigma_over_radius,11.)-pow(sigma_over_radius,5.))*(sigma_over_radius_sqrd); //for a single pair

        double total_force = m_twntyfour_epsilon*(2.0*pow(radius,-13.)-pow(radius,-7.)); //for a single pair  //-13 and -7 due to (1/r^11-1/r^5)(1/r^2) from chain rule
         //CAN ACTUALLY IMMEDIATELY CALCULATE total_force_over_r here!! just change the power!

        //ATTRACTIVE FORCE SHOULD POINT TOWARDS OTHER ATOM. REPULSIVE AWAY FROM OTHER ATOM!!!

        //find and set force components    
        double total_force_over_r = total_force/radius; //precalculate to save 2 FLOPS
        for(int j=0;j<3;j++) {
            current_atom->force[j] += total_force_over_r*displacement[j]; //i.e. Fx = (F/r)*x
            other_atom->force[j] -= current_atom->force[j]; //using Newton's 3rd law  //VERIFIED THAT EQACTLY SAME AS WRITING += -current_atom....
        }

        if(system.steps() % system.m_sample_freq ==0){
            //calculate potential energy every m_sample_freq steps
            //m_potentialEnergy += m_four_epsilon*(pow(sigma_over_radius,12.)-pow(sigma_over_radius,6));
            m_potentialEnergy += m_four_epsilon*(pow(radius,-12.)-pow(radius,-6));
        }//IF CALCULATE EVERY TIME STEP VS. EVERY 10, CPU TIME IS 60sec INSTEAD OF ~43sec
    }

    //std::cout<<"current_atom total force = " << current_atom->force <<std::endl;
 }



 //Implementing it in the above way, with only 1 calculation  per pair, made CPU time decrease from ~70sec to ~42sec for 108 atoms and 50k steps
 //1fs step size calculation

//I'VE CHECKED TO MAKE SURE THAT THE ABOVE IMPLEMENTATION GIVES SAME RESULTS AS THE BELOW OLDER IMPLEMENTATION OF CALCULATING ALL FORCES EXPLICITELY.
 //Below is previous version which loops over all atoms: so uncessesarily counts each pairwise force twice

 /*
//int i = 0;

 for(Atom *current_atom : system.atoms()) {
     //loop over the entire vector m_atoms (atoms() returns m_atoms: vector of pointers to the atoms objects

     //Note: the forces for all atoms are already reset in system.cpp calculateForces() fnc.
     //check if forces are already reset
     //std::cout << "current_atom force = " <<  current_atom->force[0] << current_atom->force[1]<< current_atom->force[2]<<std::endl;
     //std::cout << current_atom <<std::endl;    //shows that current_atom is a pointer to the location of the atom position


     //i++;
    //std::cout<<"atom number = " <<i <<std::endl;


    for(Atom *other_atom : system.atoms()){     //Calculate pairwise forces
        if(other_atom == current_atom) continue; //skip case of same atom

        //distance and vector btw objects dx should obey min. image convention: only closest distance to particle or its image is considered
        vec3 displacement;
        for(int j=0;j<3;j++){
             displacement[j] = current_atom->position[j] - other_atom->position[j];  //WHEN USE POINTERS FOR ATOMS MUST USE -> TO ACCESS THEIR PROP.'S!
             //for cases where the folded back particle will be closer than its image to a given particle
             if (displacement[j] >  system.halfsystemSize(j)) displacement[j] -= system.systemSize(j);   //systemSize(j) returns m_systemSize[j] from system class
             if (displacement[j] <= -system.halfsystemSize(j)) displacement[j] += system.systemSize(j);
         }

        //precalculate
        //declare the variables inside loop since they are only needed within the loop
        double radius = displacement.length();
        double sigma_over_radius = m_sigma/radius;
        //std::cout<<m_twntyfour_epsilon<<std::endl;
        double sigma_over_radius_sqrd = sigma_over_radius/radius;
        //double sigma_over_radius_sqrd = m_sigma/(radius*radius);
        double total_force = m_twntyfour_epsilon*(2.0*pow(sigma_over_radius,11.)-pow(sigma_over_radius,5.))*(sigma_over_radius_sqrd); //for a single pair

        //ATTRACTIVE FORCE SHOULD POINT TOWARDS OTHER ATOM. REPULSIVE AWAY FROM OTHER ATOM!!!

        //find and set force components
        double total_force_over_r = total_force/radius; //precalculate to save 2 FLOPS
        for(int j=0;j<3;j++) current_atom->force[j] += total_force_over_r*displacement[j]; //i.e. Fx = (F/r)*x

        //NOTE: THIS IS CURRENTLY DOUBLE COUNTING: NEED TO REIMPLEMENT FORCE TO SAVE 1/2 OF THE COMPUTATION'S USING NEWTON'S 3RD LAW!!
        //ALSO MAKE THE 10 THE FREQUENCY THAT THE POTENTIAL IS CALCULATED BE PASSED TO THE FUNCITON FROM MAIN.CPP SO IT CAN BE EASILY  CHANGED.
        if(system.steps() % 10 ==0){
            //calculate potential energy every 10 steps
            m_potentialEnergy += m_four_epsilon*pow(sigma_over_radius,12.)-pow(sigma_over_radius,6);
        }


    }
    //std::cout<<"current_atom total force = " << current_atom->force <<std::endl;

 }

*/

}



