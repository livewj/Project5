#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    for(Atom *atom : m_atoms) {
        if(atom->position.x() < 0){
            atom->position.setX(atom->position.x() + m_systemSize.x());
            atom->initialPosition.setX(atom->initialPosition.x() + m_systemSize.x());
            atom->positionOnNeighborlistBuild.setX(atom->positionOnNeighborlistBuild.x() + m_systemSize.x());
        }
        if(atom->position.y() < 0){
            atom->position.setY(atom->position.y() + m_systemSize.y());
            atom->initialPosition.setY(atom->initialPosition.y() + m_systemSize.y());
            atom->positionOnNeighborlistBuild.setY(atom->positionOnNeighborlistBuild.y() + m_systemSize.y());
        }
        if(atom->position.z() < 0){
            atom->position.setZ(atom->position.z() + m_systemSize.z());
            atom->initialPosition.setZ(atom->initialPosition.z() + m_systemSize.z());
            atom->positionOnNeighborlistBuild.setZ(atom->positionOnNeighborlistBuild.z() + m_systemSize.z());
        }
        if(atom->position.x() > m_systemSize.x()){
            atom->position.setX(atom->position.x() - m_systemSize.x());
            atom->initialPosition.setX(atom->initialPosition.x() - m_systemSize.x());
            atom->positionOnNeighborlistBuild.setX(atom->positionOnNeighborlistBuild.x() - m_systemSize.x());
        }
        if(atom->position.y() > m_systemSize.y()){
            atom->position.setY(atom->position.y() - m_systemSize.y());
            atom->initialPosition.setY(atom->initialPosition.y() - m_systemSize.y());
            atom->positionOnNeighborlistBuild.setY(atom->positionOnNeighborlistBuild.y() - m_systemSize.y());
        }
        if(atom->position.z() > m_systemSize.z()){
            atom->position.setZ(atom->position.z() - m_systemSize.z());
            atom->initialPosition.setZ(atom->initialPosition.z() - m_systemSize.z());
            atom->positionOnNeighborlistBuild.setZ(atom->positionOnNeighborlistBuild.z() - m_systemSize.z());
        }
    }
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
    vec3 p; //momentum
    for(Atom *atom : m_atoms)
        p += atom->velocity; //total momentum
    for(Atom *atom : m_atoms){
        atom->velocity.setX(atom->velocity.x() - p.x()/m_atoms.size()); //remove momentum equally
        atom->velocity.setY(atom->velocity.y() - p.y()/m_atoms.size());
        atom->velocity.setZ(atom->velocity.z() - p.z()/m_atoms.size());
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    for(int i = 0; i < numberOfUnitCellsEachDimension; i++)
        for(int j = 0; j < numberOfUnitCellsEachDimension; j++)
            for(int k = 0; k < numberOfUnitCellsEachDimension; k++)
                for(int l = 0; l < 4; l++)  //do this for each cell
                {
                    double x,y,z;
                    Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                    if (l == 0){
                        x = i * latticeConstant;
                        y = j * latticeConstant;
                        z = k * latticeConstant;
                    }
                    if (l == 1){
                        x = i * latticeConstant+latticeConstant/2.;
                        y = j * latticeConstant+latticeConstant/2.;
                        z = k * latticeConstant;
                    }
                    if (l == 2){
                        x = i * latticeConstant+latticeConstant/2.;
                        y = j * latticeConstant;
                        z = k * latticeConstant+latticeConstant/2.;
                    }
                    if (l == 3){
                        x = i * latticeConstant;
                        y = j * latticeConstant+latticeConstant/2.;
                        z = k * latticeConstant+latticeConstant/2.;
                    }

                    atom->position.set(x,y,z);
                    atom->initialPosition.set(x,y,z);
                    atom->resetVelocityMaxwellian(temperature);
                    m_atoms.push_back(atom);

                }

    setSystemSize(vec3(latticeConstant*numberOfUnitCellsEachDimension,
                       latticeConstant*numberOfUnitCellsEachDimension,
                       latticeConstant*numberOfUnitCellsEachDimension)); //setting the correct systemt size
/*
  //UNIFORM DISTRIBUTION OF ATOMS
    for(int i=0; i<100; i++) {
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        double x = Random::nextDouble(0, 10); // random number in the interval [0,10]
        double y = Random::nextDouble(0, 10);
        double z = Random::nextDouble(0, 10);
        atom->position.set(x,y,z);
        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);
    }
    setSystemSize(vec3(10, 10, 10)); // Remember to set the correct system size!

*/
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
