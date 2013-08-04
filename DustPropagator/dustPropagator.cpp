/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 * 
 *
 *    References
 *      
 *
 *    Notes
 *
 */

#include <fstream>
#include <limits>
#include <string>
#include <utility>

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


int main( )
{
    using namespace std;
    using namespace dust_propagator;

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



    ///////////////////////////////////////////////////////////////////////


    // Input deck

    // Values are obtained from http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html

    //Set output directory
    const string outputDirectory = "/Users/sethmichaelhirsh/Desktop/Tudat/RawData/Cases";

    //Set simulation start epoch
    const double simulationStartEpoch = 0.0;

    //Set simulation end epoch
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY*7.0;

    //Set numerical integration fixed step size in seconds
    const double fixedStepSize = 60.0 ;
   
    //Astronamical Unit in meters
    //const double AU = 1.4959787071e11;

    //Eccentricity of Earth's Orbit around Sun = 0.016710220
    const double eccentricity  = 0.1;

    //Altitude of Perigee in meters
    const double altitudeOfPerigee = 5500.0e3;

    //Set Sun equatorial radius [in m] 6.95500e8
   //for earth 6,378.137 km
   const double solarEquatorialRadius = 6378.137e3;


    const double semiMajorAxis = (solarEquatorialRadius + altitudeOfPerigee) / (1.0 - eccentricity);


//Initial state of dust particle in Keplerian coordinates
    Vector6d dustInitialStateInKeplerianElements;
   dustInitialStateInKeplerianElements( semiMajorAxisIndex ) = semiMajorAxis;
   dustInitialStateInKeplerianElements( eccentricityIndex) = eccentricity;
   dustInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians(88.0);
   dustInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians(120.0);
   dustInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians(30.0);
   dustInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians(0.0);

//Set Sun gravitational parameter [in m^3/s^2]
   //for Sun 1.32712440018e20
   //for earth 3.986004418e14
   const double sunGravitationalParameter = 3.986004418e14;


//Set J2 zonal coefficient of the Sun
   //Source http://iau-comm4.jpl.nasa.gov/EVP.pdf "J2 = (1.7 ±0.5)×10−7"
   //unnormalized version (c or s coefficients)
   const double J2 = 0.0;



////////////////////////////////////////////


//Convert dust initial state from Keplerian elements to Cartesian elements
const Vector6d dustInitialState = convertKeplerianToCartesianElements(
                dustInitialStateInKeplerianElements,sunGravitationalParameter);


/////////////////////////////////////////////

//Creates acceleration model for dust particle

//Create body object containing dust initial state and epoch
const BodyPointer dust = boost::make_shared< Body >( dustInitialState );

//Calculate gravitiation acceleration model for dust particle
//Note: the Solar radius, in addition to J2, J3, and J4 are currently set to 0
const CentralJ2J3J4GravitationalAccelerationModelPointer dustGravityModel
        = boost::make_shared<CentralJ2J3J4GravitationalAccelerationModel >(
            boost::bind( &Body::getCurrentPosition, dust),
            sunGravitationalParameter, solarEquatorialRadius, J2, 0.0, 0.0);

//Create List of acceleration models
const CartesianStateDerivativeModel6d::AccelerationModelPointerVector dustListOfAccelerations
        = boost::assign::list_of(dustGravityModel);


///////////////////////////////////////////

//Create state derivative model
const CartesianStateDerivativeModel6dPointer dustStateDerivativeModel
       = boost::make_shared< CartesianStateDerivativeModel6d >(
        dustListOfAccelerations,boost::bind( &Body::setCurrentTimeAndState, dust, _1, _2) );

//////////////////////////////////////////

//Create Runge Kutta 4 integrator and set initial conditions
RungeKutta4Integrator<double, Vector6d, Vector6d > rungeKutta4(
        boost::bind( &CartesianStateDerivativeModel6d::computeStateDerivative,
                    dustStateDerivativeModel, _1, _2),
                    0.0, (Eigen::VectorXd(6) << dustInitialState).finished() );

//Set running time to initial Epoch
double runningTime = simulationStartEpoch;

DoubleKeyTypeVectorXdValueTypeMap dustPropagationHistoryCartCoord;
DoubleKeyTypeVectorXdValueTypeMap dustPropagationHistoryKeplerCoord;
DoubleKeyTypeVectorXdValueTypeMap dustPropagationHistoryCOEKeplerCoord;

//Add initial state to propagation histories
dustPropagationHistoryCartCoord[0.0] = dustInitialState;
//dustPropagationHistoryKeplerCoord[0.0] = dustInitialStateInKeplerianElements;
dustPropagationHistoryCOEKeplerCoord[0.0] = dustInitialStateInKeplerianElements;




while (runningTime < simulationEndEpoch)
{
    //Perform integration step calling compute state derivative
    Vector6d integratedState = rungeKutta4.performIntegrationStep(fixedStepSize);

    //Get the current time
    runningTime = rungeKutta4.getCurrentIndependentVariable();

    //Add state to propagation history
    dustPropagationHistoryCartCoord[ runningTime ] = integratedState;


    // Convert propagation history into Keplerian elements
    dustPropagationHistoryKeplerCoord[runningTime] = convertCartesianToKeplerianElements(
                  integratedState,sunGravitationalParameter);


   dustPropagationHistoryCOEKeplerCoord[runningTime] = propagateKeplerOrbit(
                    dustPropagationHistoryCOEKeplerCoord[runningTime - fixedStepSize],
                    fixedStepSize,
                    sunGravitationalParameter);
}





// Record mean anomaly in propagation history rather than true anomaly

runningTime = simulationStartEpoch;
double eccentricAnomaly;
double meanAnomaly;

double eccentricAnomalyCOE;
double meanAnomalyCOE;

while (runningTime < simulationEndEpoch)
{

    dustPropagationHistoryKeplerCoord[runningTime] = convertCartesianToKeplerianElements(
                  dustPropagationHistoryCartCoord[runningTime]
                  ,sunGravitationalParameter);
   //Convert true anomaly to eccentric anomaly
       eccentricAnomaly = convertTrueAnomalyToEccentricAnomaly(
                            dustPropagationHistoryKeplerCoord[runningTime](trueAnomalyIndex),
                            dustPropagationHistoryKeplerCoord[runningTime](eccentricityIndex)
                            );

    //Convert eccentric anomaly to mean anomaly
    meanAnomaly = convertEccentricAnomalyToMeanAnomaly(
                        eccentricAnomaly,
                        dustPropagationHistoryKeplerCoord[runningTime](eccentricityIndex)
                        );


    //Record mean anomaly in history to replace to true anomaly
    dustPropagationHistoryKeplerCoord[runningTime](trueAnomalyIndex) = meanAnomaly;

    //Convert inclinationIndex, argument of pariapsis index, longitude of ascending node index, and true anomaly Index to degrees
    dustPropagationHistoryKeplerCoord[runningTime](inclinationIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryKeplerCoord[runningTime](inclinationIndex));

    dustPropagationHistoryKeplerCoord[runningTime](argumentOfPeriapsisIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryKeplerCoord[runningTime](argumentOfPeriapsisIndex));

    dustPropagationHistoryKeplerCoord[runningTime](longitudeOfAscendingNodeIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryKeplerCoord[runningTime](longitudeOfAscendingNodeIndex));

    dustPropagationHistoryKeplerCoord[runningTime](trueAnomalyIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryKeplerCoord[runningTime](trueAnomalyIndex));






//For propagation of dust particle using COE

    //Convert true anomaly to eccentric anomaly
    eccentricAnomalyCOE = convertTrueAnomalyToEccentricAnomaly(
                            dustPropagationHistoryCOEKeplerCoord[runningTime](trueAnomalyIndex),
                            dustPropagationHistoryCOEKeplerCoord[runningTime](eccentricityIndex)
                            );

    //Convert eccentric anomaly to mean anomaly
    meanAnomalyCOE = convertEccentricAnomalyToMeanAnomaly(
                        eccentricAnomalyCOE,
                        dustPropagationHistoryCOEKeplerCoord[runningTime](eccentricityIndex)
                        );


    //Record mean anomaly in history to replace to true anomaly
    dustPropagationHistoryCOEKeplerCoord[runningTime](trueAnomalyIndex) = meanAnomalyCOE;

    //Convert inclinationIndex, argument of pariapsis index, longitude of ascending node index, and true anomaly Index to degrees
    dustPropagationHistoryCOEKeplerCoord[runningTime](inclinationIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryCOEKeplerCoord[runningTime](inclinationIndex));

    dustPropagationHistoryCOEKeplerCoord[runningTime](argumentOfPeriapsisIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryCOEKeplerCoord[runningTime](argumentOfPeriapsisIndex));

    dustPropagationHistoryCOEKeplerCoord[runningTime](longitudeOfAscendingNodeIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryCOEKeplerCoord[runningTime](longitudeOfAscendingNodeIndex));

    dustPropagationHistoryCOEKeplerCoord[runningTime](trueAnomalyIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryCOEKeplerCoord[runningTime](trueAnomalyIndex)); 


    runningTime = runningTime + fixedStepSize;
}


//////////////////////////////////////

//Write DataMap to text file in Cartesian Elements
writeDataMapToTextFile( dustPropagationHistoryCartCoord, 
                        "dustPropagationHistoryDNICartCase4_60.dat",
                        outputDirectory,
                        "",
                        numeric_limits< double >::digits10,
                        numeric_limits< double >::digits10,
                        ",");

//write DataMap to text file in Keplerian Elements
writeDataMapToTextFile( dustPropagationHistoryKeplerCoord, 
                        "dustPropagationHistoryDNICase4_60s.dat",
                        outputDirectory,
                        "",
                        numeric_limits< double >::digits10,
                        numeric_limits< double >::digits10,
                        ",");


//write DataMap to text file in Keplerian Elements
writeDataMapToTextFile( dustPropagationHistoryCOEKeplerCoord, 
                        "dustPropagationHistoryCOECase4_60s.dat",
                        outputDirectory,
                        "",
                        numeric_limits< double >::digits10,
                        numeric_limits< double >::digits10,
                        ",");


    
    return 0;
}