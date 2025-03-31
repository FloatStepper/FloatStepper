/*--------------------------------------------------- -----------------------*\
|   Module Name:     FloatStepper                                             |
|   Description:     OpenFOAM extension module for fluid-rigid body coupling  |
|   License:         GNU General Public License (GPL) version 3               |
|   Copyright:       2025 Johan Roenby, STROMNING APS                         |
|---------------------------------------------------- ------------------------|
|-------Diversity-Equality-Inclusion----Slava-Ukraini----Free-Palestine-------|
\*--------------------------------------------------- -----------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::floaterMotionRestraints::linearSpring

Description
    floaterMotionRestraints model.  Linear spring.

SourceFiles
    linearSpring.C

\*---------------------------------------------------------------------------*/

#ifndef linearSpring_H
#define linearSpring_H

#include "floaterMotionRestraint.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace floaterMotionRestraints
{

/*---------------------------------------------------------------------------*\
                          Class linearSpring Declaration
\*---------------------------------------------------------------------------*/

class linearSpring
:
    public floaterMotionRestraint
{
    // Private data

        //- Anchor point, where the spring is attached to an immovable
        //  object
        point anchor_;

        //- Reference point of attachment to the solid body
        point refAttachmentPt_;

        //- Spring stiffness coefficient (N/m)
        scalar stiffness_;

        //- Damping coefficient (Ns/m)
        scalar damping_;

        //- Rest length - length of spring when no forces are applied to it
        scalar restLength_;


public:

    //- Runtime type information
    TypeName("linearSpring");


    // Constructors

        //- Construct from components
        linearSpring
        (
            const word& name,
            const dictionary& rBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<floaterMotionRestraint> clone() const
        {
            return autoPtr<floaterMotionRestraint>
            (
                new linearSpring(*this)
            );
        }


    //- Destructor
    virtual ~linearSpring();


    // Member Functions

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const floaterMotion& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& rBMRCoeff);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
