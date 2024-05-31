/*---------------------------------------------------------------------------*\
    Module Name:       FloatStepper
    Description:       OpenFOAM extension module for fluid-rigid body coupling
    License:           GNU General Public License (GPL) version 3
    Copyright:         2023 Johan Roenby, STROMNING APS
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of FloatStepper.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FloatStepper.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::floaterMotion

Description
    Six degree of freedom motion for a floating rigid body.

    Can add restraints (e.g. a spring).

SourceFiles
    floaterMotionI.H
    floaterMotion.C
    floaterMotionIO.C

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

        //- Moment of inertia of the body in reference configuration
        //  (Q = I)
        symmTensor momentOfInertia_;

        //- Added mass matrix in body reference frame
        SquareMatrix<scalar> Madd_;

        //- How often added mass should be updated
        label MaddUpdateFreq_;

        //- Force contribution from all but added mass
        vector F0_;

        //- Torque contribution from all but added mass
        vector tau0_;

        //- Restraint force
        vector restraintForce_;

        //- Restraint torque
        vector restraintTorque_;

        //- Names of patches making up the rigid body
        wordRes patches_;

    // Private Member Functions

        // Access functions retained as private because of the risk of
        // confusion over what is a body local frame vector and what is global

        // Access

            //- Return the restraints
            inline const PtrList<floaterMotionRestraint>&
                restraints() const;

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
    ~floaterMotion();


    // Member Functions

        // Access

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


        //- Write
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
