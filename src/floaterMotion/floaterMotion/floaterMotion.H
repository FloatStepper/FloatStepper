/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|---------------------------------------------------- ------------------------|
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/

#ifndef floaterMotion_H
#define floaterMotion_H

#include "floaterMotionState.H"
#include "pointField.H"
#include "floaterMotionRestraint.H"
#include "Tuple2.H"
#include "pimpleControl.H" //Used by updateBodyState()
#include "addedMass.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

/*---------------------------------------------------------------------------*\
                      Class floaterMotion Declaration
\*---------------------------------------------------------------------------*/

class floaterMotion
{

    // Private data

        // Floater name
        const word name_;

        //- Reference to time database
        const Time& time_;

        //- Motion state data object
        floaterMotionState motionState_;

        //- Motion state data object for previous time-step
        floaterMotionState motionState0_;

        //- Motion restraints
        PtrList<floaterMotionRestraint> restraints_;

        //- Motion constraints
        //PtrList<floaterMotionConstraint> constraints_;

        //- Translational constraint tensor
        //tensor tConstraints_;

        //- Rotational constraint tensor
        //tensor rConstraints_;

        //- Centre of mass of initial state
        point initialCentreOfMass_;

        //- Centre of rotation of initial state
        point initialCentreOfRotation_;

        //- Orientation of initial state
        tensor initialQ_;

        //- Mass of the body
        scalar mass_;

        //- Moment of inertia of the body with respect to centreOfRotation and
        //  represented in the body coordinate system spanned by the columns of
        //  Q. To obtain inertia tensor around CofR from inertia tensor around
        //  ICofM, use the parallel axis theorem:
        //
        //      I_CofR = I_CofM + mass_*(I*magSqr(momentArm) - sqr(momentArm));
        //
        //  To obtain a inertia tensor in body coordinates, I_body, from the 
        //  inertia tensor in lab coordinates, I_lab, employ the coordinate
        //  transformation:
        //
        //      I_body = (Q.T() & I_lab) & Q;
        symmTensor momentOfInertia_;

        // Point with respect to which the moment of inertia is given by the
        // user. Used to calculated momentOfInertia wrt centre of rotation.
        vector momentOfInertiaRefPoint_;

        // Coordinate axes relative to which the moment of inertia matrix is
        // given by the user. Used to calculated momentOfInertia wrt lab axes.
        tensor momentOfInertiaAxes_;

        //- Added mass matrix in body reference frame with respect to CofR.
        SquareMatrix<scalar> Madd_;

        //- How often added mass should be updated
        label MaddUpdateFreq_;

        //- Force contribution from all but added mass. In lab coordinates.
        vector F0_;

        //- Torque contribution from all but added  with respect to centre of
        //  rotation, CofR. In lab coordinates.
        vector tau0_;

        //- Restraint force. In lab coordinates.
        vector restraintForce_;

        //- Restraint torque with respect to centre of rotation, CofR. In lab
        //  coordinates.
        vector restraintTorque_;

        //- Names of patches making up the rigid body
        wordRes patches_;

    // Private Member Functions

        // Access functions retained as private because of the risk of
        // confusion over what is a body local frame vector and what is global

        // Access

            //- Return the restraints
            inline const PtrList<floaterMotionRestraint>& restraints() const;

            //- Return the constraints
            //inline const PtrList<floaterMotionConstraint>&
            //    constraints() const;

            //- Return the initial centre of rotation
            inline const point& initialCentreOfRotation() const;

            //- Return the initial orientation
            inline const tensor& initialQ() const;

            //- Return the orientation
            inline const tensor& Q() const;


        // Edit

            //- Return the centre of rotation
            inline point& initialCentreOfRotation();

            //- Return initial orientation
            inline tensor& initialQ();

            //- Return non-const access to the orientation
            inline tensor& Q();

public:

    //- Runtime type information
    TypeName("floaterMotion");

    // Static Data
//    static int debug;

    // Constructors

        //- Construct null
        floaterMotion(const Time&);

        //- Construct from constant and state dictionaries
        floaterMotion
        (
            const dictionary& dict,
            const dictionary& stateDict,
            const Time& time
        );

        //- Construct as copy
        floaterMotion(const floaterMotion&);


    //- Destructor
    virtual ~floaterMotion() = default;


    // Member Functions

        // Access

            //- Return the name
            inline word name() const;

            //- Return the mass
            inline scalar mass() const;

            //- Return the inertia tensor
            inline const symmTensor& momentOfInertia() const;

            //- Return the motion state
            inline const floaterMotionState& state() const;

            //- Return the current centre of rotation
            inline const point& centreOfRotation() const;

            //- Return the initial centre of mass
            inline const point& initialCentreOfMass() const;

            //- Return the current centre of mass
            inline point centreOfMass() const;

            //- Return the orientation tensor, Q.
            //  globalVector = Q & bodyLocalVector
            //  bodyLocalVector = Q.T() & globalVector
            inline const tensor& orientation() const;

            //- Return the current velocity
            inline const vector& v() const;

            //- Return the angular velocity in the global frame
            inline const vector& omega() const;

            inline vector momentArm() const;

            //- Return time
            inline const Time& time() const;

            inline const vector& restraintForce() const;

            inline const vector& restraintTorque() const;

            inline const SquareMatrix<scalar>& Madd() const;

            inline const vector& F0() const;

            inline const vector& tau0() const;


        // Edit

            //- Store the motion state at the beginning of the time-step
            inline void storeOldState();

            //- Return non-const access to the centre of rotation
            inline point& centreOfRotation();

            //- Set the current velocity JRO
            inline void setVelocity(const vector& V);

            //- Set the current angular velocity JRO
            inline void setAngularVelocity(const vector& omega);

            //- Set the current centre of rotation JRO
            inline void setCentreOfRotation(const vector& CoR);

            //- Set the current orientation JRO
            inline void setOrientation(const tensor& Q);

            //- Set the current linear acceleration JRO
            inline void setAcceleration(const vector& a);

            //- Set the current angular acceleration JRO
            inline void setAngularAcceleration(const vector& domegadt);

            //- Set the current motion state JRO
            inline void setState(const floaterMotionState& state);

            //- Calculate added mass matrix
            void calcAddedMass
            (
                const fvMesh& mesh,
                const labelList& DoFs
            );


            //- Calculate full force and torque using forceContributions 
            //  functionObject
            void calcForceAndTorque
            (
                const volScalarField& rho,
                const volScalarField& p,
                const volVectorField& U,
                vector& force,
                vector& torque
            );


           //- Integrate equations of motion JRO
           void updateFloaterState(const scalarField& dvwdt, const scalar dt);


           //- Assemble and solve body equations of motion given the added mass
           //  tensor and the force and torque contributions, F0 and tau0,
           //  representing the part of the total force and torque which is not
           //  linear in the body acceleration. 
            scalarField calcAcceleration
            (
                const vector& F0,
                const vector& tau0,
                const labelList& DoFs
            );

            //- Change frame of 6x6 matrices. For instance to go from lab to 
            //  body frame supply the body orientation matrix (with body axes in
            //  its columns) as 2nd argument. To go back again supply its
            //  transpose.
            scalarSquareMatrix changeFrame
            (
                const scalarSquareMatrix& M,
                const tensor& Q
            ) const;

            // Is it time for recalculating added mass matrix?
            bool MaddUpdateTime() const;

        // Constraints and Restraints

            //- Add restraints to the motion, public to allow external
            //  addition of restraints after construction
            void addRestraints(const dictionary& dict);

            //- Apply the restraints to the object (made public by JRO)
            void applyRestraints();

            //- Add restraints to the motion, public to allow external
            //  addition of restraints after construction
            void addConstraints(const dictionary& dict);

            //- Report the status of the motion
            void status() const;


        // Transformations

            //- Return the velocity of a position
            inline point velocity(const point& pt) const;

            //- Transform the given initial state point by the current motion
            //  state
            inline point transform(const point& initialPoints) const;

            //- Transform the given initial state pointField by the current
            //  motion state
            tmp<pointField> transform(const pointField& initialPoints) const;

            //- Transform the given initial state pointField by the current
            //  motion state scaled by the given scale
            tmp<pointField> transform
            (
                const pointField& initialPoints,
                const scalarField& scale
            ) const;


        //- Write to dictionary
        void write(dictionary& dict) const;

        //- Write to stream
        void write(Ostream&) const;

        //- Read coefficients dictionary and update system parameters,
        //  constraints and restraints but not the current state
        bool read(const dictionary& dict);


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "floaterMotionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
