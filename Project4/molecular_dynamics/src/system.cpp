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
    //use version A: where fold the atoms back into the simulation cell. The cell has its lower left corner at origin of coord. system.

    for(Atom *atom : atoms()) {
        for(int j=0;j<3;j++){
            //fold atoms back into the box if they escape the box
            //old PBCs: don't take into account that atom in 1 step can fly out >1 systemsize away from current position
            if (atom->position[j] <  0.) atom->position[j] += m_systemSize[j];  //I think use atom-->position instead of atom.position b/c atom is a pointer
            if (atom->position[j] >  m_systemSize[j]) atom->position[j] -= m_systemSize[j];

            //NOTE: WHEN I FIX THE ISSUE OF ATOMS FLYING OUT OF THE SYSTEM, I GET THE ISSUE OF VELOCITY EXPLODING!! PERHAPS THERE IS TOO MUCH ENERGY
            //SOMEHOW WHEN HAVE MORE ATOMS: TOO MANY COLLISIONS ETC...
            //I THINK NEED TO PERIODICALLY COUPLE THE SYSTEM TO THE MAXWELL BOLTZMANN DISTRIBUTION: GIVE PARTICLES NEW VELOCITIES: OTHERWISE BECAUSE OF
            //ARTIFICIALLY CONFINING THEM TO 1 SIMULATION CELL, ARE INCREASING THE ENERGY TOO MUCH!!!

            //OR JUST DON'T CONFINE THEM TO THE CELL: LET THEM EXCAPE IF GO BEYON 1 CELL AWAY WITHIN A TIME STEP!!!!

            //New PBCs: try to work for any position of atoms: NOTE SURE IF WORKING PROPERLY
            //if(atom->position[j] <  0.) atom->position[j] += m_systemSize[j]*floor(abs(atom->position[j])/m_systemSize[j]);
            //if(atom->position[j] > m_systemSize[j]) atom->position[j] -= m_systemSize[j]*floor(atom->position[j]/m_systemSize[j]);

            //ANOTHER VERSION which is probably less efficient
            //while (atom->position[j] <  0.) atom->position[j] += m_systemSize[j];  //I think use atom-->position instead of atom.position b/c atom is a pointer
            //while (atom->position[j] >  m_systemSize[j]) atom->position[j] -= m_systemSize[j];
            //while position is outside of the main simulation cell, keep moving the particle by 1 system size lenght
            //until position is back within simulation cell. This is necessary b/c sometimes particles with high KE fly
            //out of the system way past the neighboring image cells.

            /*
            //for using center of system as origin. This is  not convinient when building a lattice
            if (atom->position[j] <  -m_systemSize[j] * 0.5) atom->position[j] += m_systemSize[j];  //I think use atom-->position instead of atom.position b/c atom is a pointer
            if (atom->position[j] >=  m_systemSize[j] * 0.5) atom->position[j] -= m_systemSize[j];
            */
        }
    }
}

void System::removeTotalMomentum() {
    //temporarily insert numAtoms here: LATER PASS IT TO THE FUNCTION OR FIND IT SOMEHOW FROM SYSTEM PROP.'S
    int numAtoms = 500;
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
    vec3 total_momentum;
    for(Atom *atom : atoms()) { //c++11 way of iterating through  entire vector or array
        total_momentum += atom->mass()*atom->velocity;  //mass() returns value of m_mass (atom's mass)
    }
    vec3 Mom_per_atom;   //3D momentum components
    Mom_per_atom = total_momentum/numAtoms; //this is amount of abs(momentum) per atom that need to remove to get total to 0
    for(Atom *atom : atoms()) {
        for(int j=0;j<3;j++){
            //evenly modify components of velocity of each atom to yield total system momentum of 0
            atom->velocity[j] -= Mom_per_atom[j]/atom->mass();
            //if(atom->velocity[j] >0) atom->velocity[j] -= Mom_per_atom[j]/atom->mass();
            //if(atom->velocity[j] <0)  atom->velocity[j] += Mom_per_atom[j]/atom->mass();
        }
    }
    //test if total_momentum was rezeroed 0: (a unit test)
    total_momentum.print("Total Momentum before removeMomentum");
    total_momentum.set(0,0,0);  //reset total momentum to 0
    for(Atom *atom : atoms()) { //c++11 way of iterating through  entire vector or array
          total_momentum += atom->mass()*atom->velocity;  //mass() returns value of m_mass (atom's mass)
     }
    total_momentum.print("Total Momentum after removeMomentum");   //print() is fnc in vec3 class, takes in a string input
}



void System::createFCCLattice(vec3 numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    vec3 LatticeVector;  //vector which points to the origin of each unit cell
    //Note: 1st unit cell starts at 0,0,0

    double x,y,z;
    double halfLatticeConstant=0.5*latticeConstant;
    //std::cout<<"halflatticeconsts" << halfLatticeConstant << std::endl;

    for(int i=0;i<numberOfUnitCellsEachDimension[0];i++){
        //i.e. i = 0,1...N_x-1
        for(int j=0;j<numberOfUnitCellsEachDimension[1];j++){
            for(int k=0;k<numberOfUnitCellsEachDimension[2];k++){
                LatticeVector.set(latticeConstant*i,latticeConstant*j,latticeConstant*k);

                //Place the 4 atoms of each fcc cell into coordinates
                //NOTE: The PBCs will prevent from adding atoms which are beyond the system dimensions when approach the boundaries.
                Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26)); //uses mass in kg
                x = LatticeVector[0];
                y = LatticeVector[1];
                z = LatticeVector[2];
                atom1->position.set(x,y,z);
                atom1->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom1);     //add element to vector m_atoms 1 element (atom object)

                Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = halfLatticeConstant + LatticeVector[0];
                y = halfLatticeConstant + LatticeVector[1];
                z = LatticeVector[2];
                atom2->position.set(x,y,z);
                atom2->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom2);

                Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = LatticeVector[0];
                y = halfLatticeConstant + LatticeVector[1];
                z = halfLatticeConstant + LatticeVector[2];
                atom3->position.set(x,y,z);
                atom3->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom3);

                Atom *atom4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = halfLatticeConstant + LatticeVector[0];
                y = LatticeVector[1];
                z = halfLatticeConstant + LatticeVector[2];
                atom4->position.set(x,y,z);
                atom4->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom4);
            }
        }
    }


    setSystemSize(latticeConstant*numberOfUnitCellsEachDimension); //system size set by multiply vec3 # of unit cells by latticeConstant
    std::cout<<"system size = " << m_systemSize <<std::endl;
    //each unit cell is cubic with side length = latticeConstant

    /*
    //Places 100 atoms randomly into a cube
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
    */
}


void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt);  //this calls velocityverlet.cpp
    m_steps++;
    m_time += dt;
}
