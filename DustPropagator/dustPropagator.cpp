/*    
 *    Copyright (c) 2013, S. Hirsh
 *    Copyright (c) 2013, K. Kumar
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING or http://bit.ly/1cFRf75 for license details.
 */

#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <sstream>

#include <boost/assign/list_of.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <TudatCore/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h>
#include <Tudat/Astrodynamics/Gravitation/centralJ2J3J4GravityModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/cartesianStateDerivativeModel.h>
#include <Tudat/Astrodynamics/StateDerivativeModels/compositeStateDerivativeModel.h>
#include <Tudat/InputOutput/basicInputOutput.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

#include "DustPropagator/body.h"

#include <Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h>

// Declare Kepler elements enum with mean anomaly instead of true anomaly.
enum keplerianElementsWithMeanAnomaly
{
    semiMajorAxisIndex,
    eccentricityIndex,
    inclinationIndex,
    longitudeOfAscendingNodeIndex,
    argumentOfPeriapsisIndex,
    meanAnomalyIndex
};

int main( )
{
    using namespace std;

    using tudat::basic_astrodynamics::AccelerationModel3dPointer;
    using tudat::basic_astrodynamics::semiMajorAxisIndex;
    using tudat::basic_astrodynamics::eccentricityIndex;
    using tudat::basic_astrodynamics::inclinationIndex;
    using tudat::basic_astrodynamics::argumentOfPeriapsisIndex;
    using tudat::basic_astrodynamics::longitudeOfAscendingNodeIndex;
    using tudat::basic_astrodynamics::trueAnomalyIndex;

    using tudat::basic_mathematics::Vector6d;

    using tudat::gravitation::CentralJ2J3J4GravitationalAccelerationModel;
    using tudat::gravitation::CentralJ2J3J4GravitationalAccelerationModelPointer;

    using tudat::input_output::DoubleKeyTypeVectorXdValueTypeMap;
    using tudat::input_output::writeDataMapToTextFile;

    using tudat::numerical_integrators::RungeKutta4Integrator;

    using tudat::orbital_element_conversions::convertKeplerianToCartesianElements;
    using tudat::orbital_element_conversions::convertCartesianToKeplerianElements;

    using tudat::state_derivative_models::CartesianStateDerivativeModel6d;
    using tudat::state_derivative_models::CartesianStateDerivativeModel6dPointer;
    using tudat::state_derivative_models::CompositeStateDerivativeModel;

    using tudat::unit_conversions::convertDegreesToRadians;
    using tudat::unit_conversions::convertRadiansToDegrees;

    using tudat::orbital_element_conversions::convertTrueAnomalyToEccentricAnomaly;
    using tudat::orbital_element_conversions::convertEccentricAnomalyToMeanAnomaly;

    using tudat::orbital_element_conversions::propagateKeplerOrbit;

    using namespace dust_propagator;

    ///////////////////////////////////////////////////////////////////////

    // KK: Don't forget to comment code blocks.
    // KK: Ensure that you declare as much as possible outside the for-loops.
    // KK: Be explicit with variable-naming so the code reads easily.

    // Set up input deck.

    // Set output directory.
    const string outputDirectory = "/Users/kartikkumar/Desktop/";

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY * 7.0;

    // Set system parameters.
    // Values are obtained from http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html

    // Set central body's equatorial radius (m).
    const double centralBodyEquatorialRadius = 6378.137e3;

    // Set central body's gravitational parameter (m^3/s^2).
    // mu_Sun = 1.32712440018e20
    // mu_Earth = 3.986004418e14
    const double centralBodyGravitationalParameter = 3.986004418e14;


    // Set J2 zonal coefficient for central body.
    // Source http://iau-comm4.jpl.nasa.gov/EVP.pdf "J2 = (1.7 ±0.5)×10−7"
    // unnormalized version (c or s coefficients)
    const double J2 = 0.0;

    // 6 x 5 array containing all of the test data. Each row is a different element of the initial 
    // conditions.
    // Each column is a different case. 
    const int numberOfCases = 5;

    // KK: Most likely safer to use STL containers for this (e.g., vector of vectors).
    const double cases[ 6 ][ numberOfCases ] = 
    {
        { 360.0e3,  550.0e3, 1350.0e3, 5500.0e3, 360.0e3 }, // Altitude of Perigee (m)
        { 0.0005,   0.005,   0.05,     0.1,      0.005 },   // Eccentricity 
        { 10.0,     97.5,    101.0,    88.0,     10.0 },    // Inclination (deg)
        { 30.0,     50.0,    80.0,     120.0,    30.0 },    // Longitude of Ascending Node (deg)
        { 30.0,     30.0,    30.0,     30.0,     30.0 },    // Argument of Periapsis (deg)
        { 0.0,      0.0,     0.0,      0.0,      0.0 }      // True anomaly (deg)
    };

    // Set integrator step sizes (s).
    const int numberOfStepSizes = 2;
    const double stepSizes[ numberOfStepSizes ] = { 10.0, 60.0 };

    ///////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////

    // Execute simulations.

    // Loop over cases.
    // KK: Be very careful with array indices: in C++ array indices start at 0! Convention is to 
    // start for-loop counters at 0 and use the '<' comparison to check of the end condition.
    for ( int caseNumber = 0; caseNumber < numberOfCases; caseNumber++ )
    { 
        // Loop over numerical integrator step sizes.
        for (int stepSizeNumber = 0; stepSizeNumber < numberOfStepSizes; stepSizeNumber++ )
        {
            ///////////////////////////////////////////////////////////////////////

            // Set up initial conditions.

            // Set numerical integration fixed step size (s).
            const double fixedStepSize = stepSizes[ stepSizeNumber ];
            
            // Altitude of perigee (m).
            const double altitudeOfPerigee = cases[ 0 ][ caseNumber ];   
            
            // Eccentricity of Earth's Orbit around Sun = 0.016710220.
            const double eccentricity = cases[ 1 ][ caseNumber ];
            
            // Compute semi-major axis (m).
            const double semiMajorAxis = ( centralBodyEquatorialRadius + altitudeOfPerigee ) 
                                         / ( 1.0 - eccentricity );

            // Initial state of dust particle in Keplerian coordinates
            Vector6d dustInitialStateInKeplerianElements;
            dustInitialStateInKeplerianElements( semiMajorAxisIndex ) = semiMajorAxis;
            dustInitialStateInKeplerianElements( eccentricityIndex) = eccentricity;
            dustInitialStateInKeplerianElements( inclinationIndex ) 
                = convertDegreesToRadians( cases[ 2 ][ caseNumber ] );
            dustInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) 
                = convertDegreesToRadians( cases[ 3 ][ caseNumber ] );
            dustInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) 
                = convertDegreesToRadians( cases[ 4 ][ caseNumber ] );
            dustInitialStateInKeplerianElements( trueAnomalyIndex ) 
                = convertDegreesToRadians( cases[ 5 ][ caseNumber ] );

            // Convert dust initial state from Keplerian elements to Cartesian elements
            const Vector6d dustInitialState = convertKeplerianToCartesianElements(
                dustInitialStateInKeplerianElements, centralBodyGravitationalParameter );

            ///////////////////////////////////////////////////////////////////////

            ///////////////////////////////////////////////////////////////////////    

            // Create acceleration model(s) and state derivative model for dust particle.

            // Create body object containing dust initial state and epoch.
            const BodyPointer dust = boost::make_shared< Body >( dustInitialState );

            // Calculate gravitation acceleration model for dust particle.
            // Note: the central body radius, in addition to J2, J3, and J4 are currently set to 0.
            const CentralJ2J3J4GravitationalAccelerationModelPointer dustGravityModel
                    = boost::make_shared< CentralJ2J3J4GravitationalAccelerationModel >(
                        boost::bind( &Body::getCurrentPosition, dust ),
                        centralBodyGravitationalParameter, centralBodyEquatorialRadius, 
                        J2, 0.0, 0.0 );

            // Create list of acceleration models.
            const CartesianStateDerivativeModel6d::AccelerationModelPointerVector 
                dustListOfAccelerations = boost::assign::list_of( dustGravityModel );

            // Create state derivative model.
            const CartesianStateDerivativeModel6dPointer dustStateDerivativeModel
                   = boost::make_shared< CartesianStateDerivativeModel6d >(
                        dustListOfAccelerations,boost::bind( 
                            &Body::setCurrentTimeAndState, dust, _1, _2 ) );                

            ///////////////////////////////////////////////////////////////////////    

            ///////////////////////////////////////////////////////////////////////

            // Create Runge Kutta 4 integrator and set initial conditions.
            RungeKutta4Integrator< double, Vector6d, Vector6d > rungeKutta4(
                    boost::bind( &CartesianStateDerivativeModel6d::computeStateDerivative,
                                 dustStateDerivativeModel, _1, _2 ),
                    simulationStartEpoch, 
                    ( Eigen::VectorXd( 6 ) << dustInitialState ).finished( ) );

            /////////////////////////////////////////////////////////////////////// 

            ///////////////////////////////////////////////////////////////////////

            // Execute numerical integration.

            // Set running time to initial epoch.
            double runningTime = simulationStartEpoch;

            // Declare state history containers.
            DoubleKeyTypeVectorXdValueTypeMap dustStateHistory;
            DoubleKeyTypeVectorXdValueTypeMap dustStateHistoryKeplerElements;
            DoubleKeyTypeVectorXdValueTypeMap dustStateHistoryKeplerElementsBenchmark;

            // Add initial state to state histories.
            dustStateHistory[ runningTime ] = dustInitialState;
            dustStateHistoryKeplerElements[ runningTime ] 
                = dustInitialStateInKeplerianElements;
           dustStateHistoryKeplerElementsBenchmark[ runningTime ] 
                = dustInitialStateInKeplerianElements;

            // Integrate until the simulation end epoch.
            while (runningTime < simulationEndEpoch)
            {
                // Perform integration step calling compute state derivative.
                Vector6d integratedState = rungeKutta4.performIntegrationStep( fixedStepSize );

                // Get the current time.
                runningTime = rungeKutta4.getCurrentIndependentVariable( );

                // Add state to state history.
                dustStateHistory[ runningTime ] = integratedState;

                // Convert state history into Keplerian elements.
                dustStateHistoryKeplerElements[ runningTime ] 
                    = convertCartesianToKeplerianElements( 
                        integratedState, centralBodyGravitationalParameter );

                dustStateHistoryKeplerElementsBenchmark[ runningTime ] 
                    = propagateKeplerOrbit(
                        dustStateHistoryKeplerElementsBenchmark[ runningTime - fixedStepSize ],
                        fixedStepSize,
                        centralBodyGravitationalParameter);
            }

            ///////////////////////////////////////////////////////////////////////

            ///////////////////////////////////////////////////////////////////////

            // Convert state history data to desired output format.
            // Convert radians to degrees for orbital elements.
            // Record mean anomaly instead of true anomaly.

            // Create new data maps to record final formats of dust state histories.
            DoubleKeyTypeVectorXdValueTypeMap dustStateHistoryKeplerElementsFinalFormat;
            DoubleKeyTypeVectorXdValueTypeMap dustStateHistoryKeplerElementsBenchmarkFinalFormat;

            // KK: It's safer to loop through data maps using iterators.

            // Loop through the dust state history in Kepler elements, recording data in new data
            // maps and formatting them as specified below.
            for ( DoubleKeyTypeVectorXdValueTypeMap::iterator iteratorDustStateHistory 
                  = dustStateHistoryKeplerElements.begin( );
                  iteratorDustStateHistory != dustStateHistoryKeplerElements.end( );
                  iteratorDustStateHistory++ )
            {        
                // Create vector to store current state in final format.
                Vector6d dustCurrentStateFormatted;

                // Record semi-major axis [m]
                dustCurrentStateFormatted( semiMajorAxisIndex ) 
                    = iteratorDustStateHistory->second( semiMajorAxisIndex );

                // Record eccentricity [-]
                dustCurrentStateFormatted( eccentricityIndex ) 
                    = iteratorDustStateHistory->second( eccentricityIndex );

                // Record inclination and convert it to degrees [deg].
                dustCurrentStateFormatted( inclinationIndex ) 
                    = convertRadiansToDegrees( 
                        iteratorDustStateHistory->second( inclinationIndex ) );

                // Record longitude of ascending node and convert it to degrees [deg].
                dustCurrentStateFormatted( longitudeOfAscendingNodeIndex )
                    = convertRadiansToDegrees( 
                        iteratorDustStateHistory->second( longitudeOfAscendingNodeIndex ) );

                // Record argument of periapsis and convert it to degrees [deg].
                dustCurrentStateFormatted(argumentOfPeriapsisIndex)
                    = convertRadiansToDegrees( 
                        iteratorDustStateHistory->second( argumentOfPeriapsisIndex ) );

                // Calculate the eccentric anomaly at the current state [rad].
                const double eccentricAnomalyCurrentState
                    = convertTrueAnomalyToEccentricAnomaly(
                        iteratorDustStateHistory->second( trueAnomalyIndex ),
                        iteratorDustStateHistory->second( eccentricityIndex ) );

                // Calculate the mean anomaly at the current state [rad].
                const double meanAnomalyCurrentState
                    = convertEccentricAnomalyToMeanAnomaly(
                        eccentricAnomalyCurrentState,
                        iteratorDustStateHistory->second( eccentricityIndex ) );

                // Record the mean anomaly and convert it to degrees [deg].
                dustCurrentStateFormatted( meanAnomalyIndex )
                    = convertRadiansToDegrees( meanAnomalyCurrentState );

                // Store the vector containing the direct 6 orbital rlements in state history
                dustStateHistoryKeplerElementsFinalFormat[ iteratorDustStateHistory->first ] 
                    = dustCurrentStateFormatted;
            }

            // Loop through the dust state history in Kepler elements for the benchmark, recording 
            // data in new data maps and formatting them as specified below.
            for ( DoubleKeyTypeVectorXdValueTypeMap::iterator iteratorDustStateHistory 
                  = dustStateHistoryKeplerElementsBenchmark.begin( );
                  iteratorDustStateHistory != dustStateHistoryKeplerElementsBenchmark.end( );
                  iteratorDustStateHistory++ )
            {        
                // Create vector to store current benchmark state in final format.
                Vector6d dustCurrentStateFormattedBenchmark;

                // Record semi-major axis [m].
                dustCurrentStateFormattedBenchmark( semiMajorAxisIndex ) 
                = iteratorDustStateHistory->second( semiMajorAxisIndex );

                // Record eccentricity [m].
                dustCurrentStateFormattedBenchmark( eccentricityIndex ) 
                    = iteratorDustStateHistory->second( eccentricityIndex );

                // Record inclination and convert it to degrees [deg].
                dustCurrentStateFormattedBenchmark( inclinationIndex ) 
                    = convertRadiansToDegrees(
                        iteratorDustStateHistory->second( inclinationIndex ) );

                // Record longitude of ascending node and convert it to degrees [deg].
                dustCurrentStateFormattedBenchmark( longitudeOfAscendingNodeIndex )
                    = convertRadiansToDegrees( 
                        iteratorDustStateHistory->second( longitudeOfAscendingNodeIndex ) );

                // Record argument of periapsis and convert it to degrees [deg].
                dustCurrentStateFormattedBenchmark( argumentOfPeriapsisIndex )
                    = convertRadiansToDegrees(
                        iteratorDustStateHistory->second( argumentOfPeriapsisIndex ) );

                // Calculate the eccentric anomaly at the current state [rad].
                const double eccentricAnomalyCurrentStateBenchmark
                    = convertTrueAnomalyToEccentricAnomaly(
                        iteratorDustStateHistory->second( trueAnomalyIndex ),
                        iteratorDustStateHistory->second( eccentricityIndex ) );

                // Calculate the mean anomaly at the current state [rad].
                const double meanAnomalyCurrentStateBenchmark
                    = convertEccentricAnomalyToMeanAnomaly(
                        eccentricAnomalyCurrentStateBenchmark,
                        iteratorDustStateHistory->second( eccentricityIndex ) );

                // Record the mean anomaly and convert it to degrees [deg].
                dustCurrentStateFormattedBenchmark( meanAnomalyIndex )
                    = convertRadiansToDegrees( meanAnomalyCurrentStateBenchmark );

                // Store the vector containing the direct 6 orbital rlements in state history.
                dustStateHistoryKeplerElementsBenchmarkFinalFormat[ 
                    iteratorDustStateHistory->first ] 
                        = dustCurrentStateFormattedBenchmark;
            }

            ///////////////////////////////////////////////////////////////////////

            ///////////////////////////////////////////////////////////////////////////////////

            // Write state histories to files.

            // Convert case number and stepsize values as strings.
            std::ostringstream fileNameBase;
            fileNameBase << "_case" << caseNumber + 1 << "_" << fixedStepSize << "s.dat";

            // Write dust state history in Cartesian elements to file.
            writeDataMapToTextFile( dustStateHistory, 
                                    "dustStateHistory" + fileNameBase.str( ),
                                    outputDirectory,
                                    "",
                                    numeric_limits< double >::digits10,
                                    numeric_limits< double >::digits10,
                                    "," );

            // Write dust state history in Keplerian elements to file.
            writeDataMapToTextFile( dustStateHistoryKeplerElements, 
                                    "dustStateHistoryKeplerElements" + fileNameBase.str( ),
                                    outputDirectory,
                                    "",
                                    numeric_limits< double >::digits10,
                                    numeric_limits< double >::digits10,
                                    ",");


            // Write dust state history in Keplerian elements for benchmark to file.
            writeDataMapToTextFile( dustStateHistoryKeplerElements, 
                                    "dustStateHistoryKeplerElementsBenchmark" 
                                    + fileNameBase.str( ),
                                    outputDirectory,
                                    "",
                                    numeric_limits< double >::digits10,
                                    numeric_limits< double >::digits10,
                                    ",");

            /////////////////////////////////////////////////////////////////////// 

        } // for-loop over step sizes
    } // for-loop over cases

    ///////////////////////////////////////////////////////////////////////
    
    // If this point is reached, exit program with success integer.
    return EXIT_SUCCESS;
}
