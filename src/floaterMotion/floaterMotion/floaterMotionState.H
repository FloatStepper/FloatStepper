/*---------------------------------------------------------------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|---------------------------------------------------- ------------------------|
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-----------------------------------------------------------------------------
License
    This file is part of FloatStepper.

    FloatStepper is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FloatStepper is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FloatStepper.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::floaterMotionState

Description
    Holds the motion state of sixDoF object.  Wrapped up together
    to allow rapid scatter to other processors.  The processors must all
    maintain exactly the same state data to avoid any drift or inconsistency.

SourceFiles
    floaterMotionStateI.H
    floaterMotionState.C
    floaterMotionStateIO.C

\*---------------------------------------------------------------------------*/

#ifndef floaterMotionState_H
#define floaterMotionState_H

#include "vector.H"
#include "point.H"
#include "symmTensor.H"
#include "tensor.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;
class floaterMotionState;
Istream& operator>>(Istream&, floaterMotionState&);
Ostream& operator<<(Ostream&, const floaterMotionState&);


/*---------------------------------------------------------------------------*\
                  Class floaterMotionState Declaration
\*---------------------------------------------------------------------------*/

class floaterMotionState
{
    // Private Data

        //- Current position of the centre of mass of the body
        point centreOfRotation_;

        //- Orientation, stored as the rotation tensor to transform
        //- from the body to the global reference frame
        //  i.e.:
        //  globalVector = Q_ & bodyLocalVector
        //  bodyLocalVector = Q_.T() & globalVector
        tensor Q_;

        //- Linear velocity of body
        vector v_;

        //- Linear acceleration of body
        vector a_;

        //- Angular velocity
        vector omega_;

        //- Angular acceleration
        vector domegadt_;


public:

    // Constructors

        //- Default construct, zero-initialized with identity transformation
        floaterMotionState();

        //- Construct from dictionary
        floaterMotionState(const dictionary& dict);

        //- Construct as copy
        floaterMotionState(const floaterMotionState&);

        //- Explicitly default the copy assignment operator
        floaterMotionState& operator=(const floaterMotionState&) = default;

    //- Destructor
    ~floaterMotionState();


    // Member Functions

        // Access

            //- Return access to the centre of mass
            inline const point& centreOfRotation() const;

            //- Return access to the orientation
            inline const tensor& Q() const;

            //- Return access to velocity
            inline const vector& v() const;

            //- Return access to acceleration
            inline const vector& a() const;

            //- Return access to angular velocity
            inline const vector& omega() const;

            //- Return access to angular acceleration
            inline const vector& domegadt() const;


        // Edit

            //- Return non-const access to the centre of mass
            inline point& centreOfRotation();

            //- Return non-const access to the orientation
            inline tensor& Q();

            //- Return non-const access to velocity
            inline vector& v();

            //- Return non-const access to acceleration
            inline vector& a();

            //- Return non-const access to angular velocity
            inline vector& omega();

            //- Return non-const access to angular acceleration
            inline vector& domegadt();


        //- Write to dictionary
        void write(dictionary& dict) const;

        //- Write to stream
        void write(Ostream& os) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, floaterMotionState&);
        friend Ostream& operator<<(Ostream&, const floaterMotionState&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "floaterMotionStateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
