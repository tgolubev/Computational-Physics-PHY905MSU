#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System()   //constructor
{

}

System::~System()  //desctructor: remove all the atoms
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();    //vector containing all atoms objects
}

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    //use version A: where fold the atoms back into the simulation cell. The cell is centered about the origin.
    //applies PBC in all directions

    for(Atom *atom : atoms()) {
        for(int j=0;j<3;j++){
            //fold atoms back into the box if they escape the box
            if (atom->position[j] <  -m_systemSize[j] * 0.5) atom->position[j] += m_systemSize[j];
            if (atom->position[j] >=  m_systemSize[j] * 0.5) atom->position[j] -= m_systemSize[j];
        }
    }



    //distance and vector btw objects dx should obey min. image convention
    //dx = x(j) - x(i)
    //if (dx >  x_size * 0.5) dx = dx - x_size;
    //if (dx <= -x_size * 0.5) dx = dx + x_size;

}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).

    for(int i=0; i<100; i++) {      //for all atoms (100 right now)
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26)); //uses mass in kg
        //Choose random x,y,z positions
        double x = Random::nextDouble(0, 10); // random number in the interval [0,10]
        double y = Random::nextDouble(0, 10);
        double z = Random::nextDouble(0, 10);
        atom->position.set(x,y,z);    //set atom to position x,y,z
        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);     //add element to vector m_atoms 1 element (atom object)
    }
    setSystemSize(vec3(10, 10, 10)); // Remember to set the correct system size! Dimensions of simulation cell. This fnc. sets m_systemSize=this vec3 that's passed
}

void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt);
    m_steps++;
    m_time += dt;
}
