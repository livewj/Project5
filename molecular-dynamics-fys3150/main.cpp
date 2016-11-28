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
#include <cmath>

using namespace std;

int main(int numberOfArguments, char **argumentList)
{

    int numberOfUnitCells = 5; //number of unit cells in each dimension
    double initialTemperature = UnitConverter::temperatureFromSI(300.0); // measured in Kelvin
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // b, measured in Angstrom


    // If a first argument is provided, it is the number of unit cells
    if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);
    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    if(numberOfArguments > 3) latticeConstant = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));

    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds.

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;

    System system;
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

    for(int timestep = 0; timestep < 1000; timestep ++) {
        system.step(dt);
        statisticsSampler.sample(system);  //sample to file statistics.txt
        if( timestep % 100 == 0 ) {
            // Print the timestep every 100 timesteps
            cout << setw(20) << system.steps() <<
                    setw(20) << system.time() <<
                    setw(20) << UnitConverter::temperatureToSI(statisticsSampler.temperature()) <<
                    setw(20) << statisticsSampler.kineticEnergy() <<
                    setw(20) << statisticsSampler.potentialEnergy() <<
                    setw(20) << statisticsSampler.totalEnergy() <<
                    setw(20) << UnitConverter::diffusionToSI( statisticsSampler.diffusionconstant() ) << endl;
        }
        movie.saveState(system);
    }

    movie.close();

    return 0;
}

/*
    int numberOfUnitCells = 5; //number of unit cells in each dimension
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // b, measured in Angstroms
    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds.

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;


    //LOOP OVER INITIAL TEMPS
    double start_init = UnitConverter::temperatureFromSI( 1700 ); //measured in Kelvin
    double max_init = UnitConverter::temperatureFromSI( 25000 );

    for (double initialTemperature = start_init; initialTemperature < max_init; initialTemperature += UnitConverter::temperatureFromSI( 50 )) {
        System system;
        system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature);
        system.potential().setEpsilon(1.0);
        system.potential().setSigma(3.405); //Angstrom

        system.removeTotalMomentum();

        StatisticsSampler statisticsSampler;

        double T_sum = 0.0;
        double T_mean;
        double D = 0;
        double N = 3000.0;

        for(int timestep=0.0; timestep < N; timestep ++ ) {
            system.step(dt);
            statisticsSampler.sample(system);  //sample to file statistics.txt
            D += statisticsSampler.diffusionconstant();
            if (timestep > 300) {  //temperature stabilized after 250 timesteps
                T_sum += statisticsSampler.temperature();
            }
        }
        T_mean = T_sum/2700.0; //Final temperature after reached equilibrium
        cout << UnitConverter::temperatureToSI( T_mean) << "   "<< UnitConverter::diffusionToSI( D/N ) << endl;

    }
    return 0;
}

*/
