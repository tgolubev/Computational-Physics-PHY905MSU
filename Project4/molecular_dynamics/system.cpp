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

    //I THINK NEED TO FOLD BACK ATOMS 1ST BEFORE CHECK MIN. IMAGE CONVENTION. BUT IN THAT CASE ARE CONFINING ALL ATOMS TO THE BOX,
    //SO HOW DO I FIND THE IMAGES!!

    //THESE PBCS ARE BASED ON USING THE CENTER AS ORIGIN: THIS IS NOT CONVENIENT: REDO TO USE LEFT CORNER!

    for(Atom *atom : atoms()) {
        for(int j=0;j<3;j++){
            //fold atoms back into the box if they escape the box
            if (atom->position[j] <  0.) atom->position[j] += m_systemSize[j];  //I think use atom-->position instead of atom.position b/c atom is a pointer
            if (atom->position[j] >  m_systemSize[j]) atom->position[j] -= m_systemSize[j];
            /*
            //for using center of system as origin. This is  not convinient when building a lattice
            if (atom->position[j] <  -m_systemSize[j] * 0.5) atom->position[j] += m_systemSize[j];  //I think use atom-->position instead of atom.position b/c atom is a pointer
            if (atom->position[j] >=  m_systemSize[j] * 0.5) atom->position[j] -= m_systemSize[j];
            */
        }
    }
    /*
    //distance and vector btw objects dx should obey min. image convention: only closest distance to particle or its image is considered
    for(Atom *atom:atoms()){
        Atom *current_atom = m_atoms[0]; //use pointer to current_atom
        for(Atom *atom: atoms()){ //THIS MUST BE CHANGED TO NOT DOUBLE COUNT: GO TO 1/2 OF ATOMS
            Atom *other_atom = m_atoms[0];
            vec3 distance;
            for(int j=0;j<3;j++){
                distance[j] = other_atom->position[j] - current_atom->position[j];  //WHEN USE POINTERS FOR ATOMS MUST USE -> TO ACCESS THEIR PROP.'S!
                //for cases where the folded back particle will be closer than its image to a given particle
                if (distance[j] >  m_systemSize[j] * 0.5) distance[j] -= m_systemSize[j];
                if (distance[j] <= -m_systemSize[j] * 0.5) distance[j] += m_systemSize[j];
            }
         }
    }
    */

}

void System::removeTotalMomentum() {
    //temporarily insert numAtoms here: LATER PASS IT TO THE FUNCTION OR FIND IT SOMEHOW FROM SYSTEM PROP.'S
    int numAtoms = 100;
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
            if(atom->velocity[j] >0) atom->velocity[j] -= Mom_per_atom[j]/atom->mass();
            if(atom->velocity[j] <0)  atom->velocity[j] += Mom_per_atom[j]/atom->mass();
        }
    }
}



void System::createFCCLattice(vec3 numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    vec3 LatticeVector;  //vector which points to the origin of each unit cell
    //Note: 1st unit cell starts at 0,0,0

    double x,y,z;
    double halfLatticeConstant=0.5*latticeConstant;
    //std::cout<<"halflatticeconsts" << halfLatticeConstant << std::endl;

    for(int i=0;i<=numberOfUnitCellsEachDimension[0];i++){
        //i.e. i = 0,1...N_x-1
        for(int j=0;j<=numberOfUnitCellsEachDimension[1];j++){
            for(int k=0;k<=numberOfUnitCellsEachDimension[2];k++){
                LatticeVector.set(latticeConstant*i,latticeConstant*j,latticeConstant*k);
                //NOTE: the resetVelocityMaxwellian currently causes the numbers to differ from set numbers, i.e. 0's become e-5
                //haven't yet figured out exactly what this is doing

                //Place the 4 atoms of each fcc cell into coordinates
                //NOTE: The PBCs will prevent from adding atoms which are beyond the system dimensions when approach the boundaries.
                Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26)); //uses mass in kg
                x = LatticeVector[0];
                y = LatticeVector[1];
                z = LatticeVector[2];
                //std::cout<<"Atom1 position = " <<x << y<< z << std::endl;
                atom1->position.set(x,y,z);
                //atom1->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom1);     //add element to vector m_atoms 1 element (atom object)

                if(i!=numberOfUnitCellsEachDimension[0] && j!=numberOfUnitCellsEachDimension[1]){
                    //These if statements are to not add atoms past the system boundaries.
                    //Note: If were to just use PBCs to loop this added atoms around, this results in
                    //having 2 atoms at 1 position at the left/bottom system boundaries.
                    Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                    x = halfLatticeConstant + LatticeVector[0];
                    y = halfLatticeConstant + LatticeVector[1];
                    z = LatticeVector[2];
                    //std::cout<<"Atom2 position = " <<x << y<< z << std::endl;
                    atom2->position.set(x,y,z);
                    //atom2->resetVelocityMaxwellian(temperature);
                    m_atoms.push_back(atom2);
                    //std::cout<<"Atom2 position = " <<atom2->position[0] <<atom2->position[1]<<atom2->position[2] << std::endl;
                    }

                if(j!=numberOfUnitCellsEachDimension[1] && k!=numberOfUnitCellsEachDimension[2]){
                    Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                    x = LatticeVector[0];
                    y = halfLatticeConstant + LatticeVector[1];
                    z = halfLatticeConstant + LatticeVector[2];
                    //std::cout<<"Atom3 position = " <<x << y<< z << std::endl;
                    atom3->position.set(x,y,z);
                    //atom3->resetVelocityMaxwellian(temperature);
                    m_atoms.push_back(atom3);
                }

                if(i!=numberOfUnitCellsEachDimension[0] && k!=numberOfUnitCellsEachDimension[2]){
                    Atom *atom4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                    x = halfLatticeConstant + LatticeVector[0];
                    y = LatticeVector[1];
                    z = halfLatticeConstant + LatticeVector[2];
                    //std::cout<<"Atom4 position = " <<x << y<< z << std::endl;
                    atom4->position.set(x,y,z);
                    // atom4->resetVelocityMaxwellian(temperature);
                    m_atoms.push_back(atom4);
                }


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
