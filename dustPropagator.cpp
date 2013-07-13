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


    ///////////////////////////////////////////////////////////////////////


    // Input deck

    // Values are obtained from http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html

    //Set output directory
    const string outputDirectory = "/Users/sethmichaelhirsh/Desktop/Tudat";

    //Set simulation start epoch
    const double simulationStartEpoch = 0.0;

    //Set simulation end epoch
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_YEAR*100.0;

    //Set numerical integration fixed step size
    const double fixedStepSize = tudat::physical_constants::JULIAN_DAY * 0.25 ;
   
    //Astronamical Unit in meters
    const double AU = 1.4959787071e11;

    //Eccentricity of Earth's Orbit
    const double eccentricity  = 0.016710220;


//Initial state of dust particle in Keplerian coordinates
    Vector6d dustInitialStateInKeplerianElements;
   dustInitialStateInKeplerianElements( semiMajorAxisIndex ) = AU;
   dustInitialStateInKeplerianElements( eccentricityIndex) = eccentricity;
   dustInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians(30.0);
   dustInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians(10.0);
   dustInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians(20.0);
   dustInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians(0.0);

//Set Sun gravitational parameter [in m^3/s^2]
   const double sunGravitationalParameter = 1.32712440018e20;

//Set Sun equatorial radius [in m] 6.95500e8
   const double solarEquatorialRadius = 6.95500e8;

//Set J2 zonal coefficient of the Sun
   //Source http://iau-comm4.jpl.nasa.gov/EVP.pdf "J2 = (1.7 ±0.5)×10−7"
   //unnormalized version (c or s coefficients)
   const double J2 = 1.7e-7;



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

DoubleKeyTypeVectorXdValueTypeMap dustPropagationHistory;
DoubleKeyTypeVectorXdValueTypeMap dustPropagationHistoryInKeplerianCoordinates;

//Add initial state to propagation histories
dustPropagationHistory[0.0] = dustInitialState;
dustPropagationHistoryInKeplerianCoordinates[0.0] = convertCartesianToKeplerianElements(
            dustInitialState,sunGravitationalParameter);

//Convert inclinationIndex, argument of pariapsis index, longitude of ascending node dndex, and true anomaly Index to degrees
    dustPropagationHistoryInKeplerianCoordinates[0.0](inclinationIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryInKeplerianCoordinates[runningTime](inclinationIndex));

    dustPropagationHistoryInKeplerianCoordinates[0.0](argumentOfPeriapsisIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryInKeplerianCoordinates[runningTime](argumentOfPeriapsisIndex));

    dustPropagationHistoryInKeplerianCoordinates[0.0](longitudeOfAscendingNodeIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryInKeplerianCoordinates[runningTime](longitudeOfAscendingNodeIndex));

    dustPropagationHistoryInKeplerianCoordinates[0.0](trueAnomalyIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryInKeplerianCoordinates[runningTime](trueAnomalyIndex));    



while (runningTime < simulationEndEpoch)
{
    //Perform integration step calling compute state derivative
    Vector6d integratedState = rungeKutta4.performIntegrationStep(fixedStepSize);

    //Get the current time
    runningTime = rungeKutta4.getCurrentIndependentVariable();

    //Add state to propagation history
    dustPropagationHistory[ runningTime ] = integratedState;


    // Convert propagation history into Keplerian elements
    dustPropagationHistoryInKeplerianCoordinates[runningTime] = convertCartesianToKeplerianElements(
                    dustPropagationHistory[runningTime],sunGravitationalParameter);

    //Convert inclinationIndex, argument of pariapsis index, longitude of ascending node dndex, and true anomaly Index to degrees
    dustPropagationHistoryInKeplerianCoordinates[runningTime](inclinationIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryInKeplerianCoordinates[runningTime](inclinationIndex));

    dustPropagationHistoryInKeplerianCoordinates[runningTime](argumentOfPeriapsisIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryInKeplerianCoordinates[runningTime](argumentOfPeriapsisIndex));

    dustPropagationHistoryInKeplerianCoordinates[runningTime](longitudeOfAscendingNodeIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryInKeplerianCoordinates[runningTime](longitudeOfAscendingNodeIndex));

    dustPropagationHistoryInKeplerianCoordinates[runningTime](trueAnomalyIndex) = convertRadiansToDegrees(
                        dustPropagationHistoryInKeplerianCoordinates[runningTime](trueAnomalyIndex));    

}


//////////////////////////////////////

//Write DataMap to text file in Cartesian Elements
writeDataMapToTextFile( dustPropagationHistory, 
                        "dustPropagationHistoryWithJ2_100Years.dat",
                        outputDirectory,
                        "",
                        numeric_limits< double >::digits10,
                        numeric_limits< double >::digits10,
                        ",");

//write DataMap to text file in Keplerian Elements
writeDataMapToTextFile( dustPropagationHistoryInKeplerianCoordinates, 
                        "dustPropagationHistoryInKeplerianCoordinatesWithJ2_100Years.dat",
                        outputDirectory,
                        "",
                        numeric_limits< double >::digits10,
                        numeric_limits< double >::digits10,
                        ",");


    
    return 0;
}
