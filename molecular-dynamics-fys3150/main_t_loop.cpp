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

    for (double initialTemperature = start_init; initialTemperature < max_init; initialTemperature += UnitConverter::temperatureFromSI( 50 )) {
        System system;
        system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature);
        system.potential().setEpsilon(1.0);
        system.potential().setSigma(3.405); //Angstrom

        system.removeTotalMomentum();

        StatisticsSampler statisticsSampler;

        double T_sum = 0.0;
        double T_mean;
        double numberOfTimesteps = 1000.0;

        for(double timestep=0.0; timestep<numberOfTimesteps; timestep ++ ) {
            system.step(dt);
            statisticsSampler.sample(system);  //sample to file statistics.txt
            if (timestep > 300 ); {
                T_sum += statisticsSampler.temperature(); //Reached stable temperature (equilibrium)
            }

        }
        T_mean = T_sum/700.0;
        cout << T_mean/initialTemperature << endl; //Final T for every T_init

    }



    return 0;
}