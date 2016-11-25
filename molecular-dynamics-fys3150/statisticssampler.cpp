#include "system.h"
#include "statisticssampler.h"
#include "lennardjones.h"
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include "unitconverter.h"
using namespace std;
//using std::ofstream; using std::cout; using std::endl;

StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::saveToFile(System &system, string filename)
{
    // Save the statistical properties for each timestep for plotting etc.
    // First, open the file if it's not open already    
    if(!m_file.good()) {
        m_file.open(filename.c_str(), ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }

    // Write to file for each timestep here
        m_file << setprecision(15) <<  system.steps() << "      "
               << setprecision(15) << system.time() << "      "
               << setprecision(15) << UnitConverter::temperatureToSI( temperature()) << "      "
               << setprecision(15) << kineticEnergy() << "      "
               << setprecision(15) << potentialEnergy() << "      "
               << setprecision(15) << totalEnergy() << "       "
               << setprecision(15) << UnitConverter::diffusionToSI( diffusionconstant() )<< endl;
}

void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    sampleDiffusionConstant(system);

    saveToFile(system, "statistics.txt"); //FILENAME HERE
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
    double v2 = 0;       //velocity squared
    for(Atom *atom : system.atoms()) {
        v2 += atom->velocity.lengthSquared();  //sum all the velocities
    }
    Atom* a = system.atoms()[1]; //extracting the total mass
    m_kineticEnergy = 0.5*a->mass()*v2;
}

void StatisticsSampler::samplePotentialEnergy(System &system)
{
    m_potentialEnergy = system.potential().potentialEnergy();
}

void StatisticsSampler::sampleTemperature(System &system)
{
    // Hint: reuse the kinetic energy that we already calculated
    m_temperature = (2.0/3)*m_kineticEnergy/(system.atoms().size());
}

void StatisticsSampler::sampleDensity(System &system)
{

}

void StatisticsSampler::sampleDiffusionConstant(System &system)
{
    double r2 = 0.0;
    m_DiffusionConstant = 0;
    //calculate the mean square displacement
    for (Atom *atom : system.atoms()){
        r2 += (atom->position.x() - atom->initialPosition.x())*(atom->position.x() - atom->initialPosition.x())
                + (atom->position.y() - atom->initialPosition.y())*(atom->position.y() - atom->initialPosition.y())
                + (atom->position.z() - atom->initialPosition.z())*(atom->position.z() - atom->initialPosition.z());
    }
    m_DiffusionConstant = r2/(6.0*system.atoms().size()*system.time());
}


