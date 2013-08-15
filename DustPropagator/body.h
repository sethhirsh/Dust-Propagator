/*    
 *    Copyright (c) 2013, S. Hirsh
 *    Copyright (c) 2013, K. Kumar
 *    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *    See COPYING or http://bit.ly/1cFRf75 for license details.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      121030    K. Kumar          File created.
 *      130225    K. Kumar          Updated include-guard and namespace names; updated Vector6d
 *                                  references to use Tudat definition.
 */

#ifndef DUST_PROPAGATOR_BODY_H
#define DUST_PROPAGATOR_BODY_H

#include <map>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

namespace dust_propagator
{

//! Test body class.
/*!
 * This class serves as an example of how a container can be constructed that stores state and
 * time information, which can be used in conjunction with acceleration models, the Cartesian state
 * derivative model, and the composite state derivative model. It should be noted that this class
 * should NOT be used "as is", without consideration for the application at hand. Classes such as
 * this are application-specific, hence unavailable through the Tudat libraries.
 */
class Body
{
public:

    //! Constructor taking a state and a time.
    /*!
     * Constructor taking an input state and time. The input state is used internally to
     * set the current position (taken as a segment of the input state given by the indices
     * (0, 3)) and the current velocity (taken as a segment of the input state given by the indices
     * (3, 3).
     * \param aState An input state vector.
     * \param aTime An input time (default = 0.0) [s].
     */
    Body( const tudat::basic_mathematics::Vector6d& aState, const double aTime = 0.0 )
        : currentState( aState ),
          currentPosition( aState.segment( 0, 3 ) ),
          currentVelocity( aState.segment( 3, 3 ) ),
          currentTime( aTime )
    { }

    //! Set current time and state.
    /*!
     * Sets the current time, position and current velocity internally based on the input
     * arguments. The current position is taken as a segment of the input state given by the
     * indices (0, 3)), and the current velocity is taken as a segment of the input state given by
     * the indices (3, 3).
     * \param aTime An input time [s].
     * \param aState An input state vector.
     */
    void setCurrentTimeAndState( const double aTime,
                                 const tudat::basic_mathematics::Vector6d& aState )
    {
        currentTime = aTime;
        currentState = aState;
        currentPosition = aState.segment( 0, 3 );
        currentVelocity = aState.segment( 3, 3 );
    }

    //! Get current state.
    /*!
     * Returns the internally stored current state vector.
     * \return Current state.
     */
    tudat::basic_mathematics::Vector6d getCurrentState( ) { return currentState; }

    //! Get current position.
    /*!
     * Returns the internally stored current position vector.
     * \return Current position.
     */
    Eigen::Vector3d getCurrentPosition( ) { return currentPosition; }

    //! Get current velocity.
    /*!
     * Returns the internally stored current velocity vector.
     * \return Current velocity.
     */
    Eigen::Vector3d getCurrentVelocity( ) { return currentVelocity; }

    //! Get current time.
    /*!
     * Returns the internally stored current time.
     * \return Current time.
     */
    double getCurrentTime( ) { return currentTime; }

protected:

private:

    //! Current state.
    tudat::basic_mathematics::Vector6d currentState;

    //! Current position.
    Eigen::Vector3d currentPosition;

    //! Current position.
    Eigen::Vector3d currentVelocity;

    //! Current time.
    double currentTime;
};

//! Typedef for shared-pointer to body.
typedef boost::shared_ptr< Body > BodyPointer; 


} // namespace dust_propagator

#endif // DUST_PROPAGATOR_BODY_H
