#include "velocityverlet.h"
#include "system.h"
#include "atom.h"

void VelocityVerlet::integrate(System &system, double dt)
{

    double half_dt = 0.5*dt;  //VERIFIED THAT GIVES EXACTLY EQUIVALENT RESULTS TO WRITING 0.5* INSIDE THE LOOPS
    if(m_firstStep) {
        system.calculateForces();
        m_firstStep = false;
    }

    for(Atom *atom : system.atoms()) {
        //this operates on the vectors directly using vec3 class
        atom->velocity += atom->force*half_dt/atom->mass();
        atom->position += atom->velocity*dt;  //NOTE: since v is computed 1st, this is actually v = vt+0.5at^2 since v = v+0.5at
    }

    system.applyPeriodicBoundaryConditions();
    system.calculateForces(); // New positions, recompute forces

    for(Atom *atom : system.atoms()) {
        atom->velocity += atom->force*half_dt/atom->mass();  //calculate new velocities
    }
}
