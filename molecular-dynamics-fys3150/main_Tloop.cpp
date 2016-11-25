#include "math/random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main(int numberOfArguments, char **argumentList)
{
    int numberOfUnitCells = 5; //number of unit cells in each dimension
    //double initialTemperature = UnitConverter::temperatureFromSI(300.0); // measured in Kelvin
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // b, measured in Angstroms
    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds.

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;



    //LOOP OVER INITIAL TEMPS
    double start_init = UnitConverter::temperatureFromSI( 1 ); //measured in Kelvin
    double max_init = UnitConverter::temperatureFromSI( 1000 );

    for (int initialTemperature = start_init; initialTemperature < max_init; initialTemperature ++) {        System system;
        system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature);
        system.potential().setEpsilon(1.0);
        system.potential().setSigma(3.405); //Angstrom

        system.removeTotalMomentum();

        StatisticsSampler statisticsSampler;

        IO movie("movie.xyz"); // To write the state to file

        cout << setw(20) << "Timestep" <<
                setw(20) << "Time" <<
                setw(20) << "Temperature" <<
                setw(20) << "KineticEnergy" <<
                setw(20) << "PotentialEnergy" <<
                setw(20) << "TotalEnergy" <<
                setw(20) << "DiffusionConstant" << endl;

        double T_sum = 0;
        double T_mean;
        double numberOfTimesteps = 1000;

        for(int timestep=0; timestep<numberOfTimesteps; timestep++) {
            system.step(dt);
            statisticsSampler.sample(system);  //sample to file statistics.txt
            T_sum += statisticsSampler.temperature();
            if( timestep % 100 == 0 ) {
                // Print the timestep every 100 timesteps
                cout << setw(20) << system.steps() <<
                        setw(20) << system.time() <<
                        setw(20) << UnitConverter::temperatureToSI(statisticsSampler.temperature()) << //convert back to SI
                        setw(20) << statisticsSampler.kineticEnergy() <<
                        setw(20) << statisticsSampler.potentialEnergy() <<
                        setw(20) << statisticsSampler.totalEnergy() <<
                        setw(20) << statisticsSampler.diffusionconstant() << endl;
            }
            movie.saveState(system);
            T_mean = T_sum/numberOfTimesteps;
            cout << UnitConverter::temperatureToSI( T_mean ) << endl; //Final T for every T_init
        }

    }

    movie.close();

    return 0;
}
