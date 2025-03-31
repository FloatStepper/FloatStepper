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
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Namespace
    Foam::floaterMotionRestraints

Description
    Namespace for six DoF motion restraints


Class
    Foam::floaterMotionRestraint

Description
    Base class for defining restraints for sixDoF motions

SourceFiles
    floaterMotionRestraint.C
    floaterMotionRestraintNew.C

\*---------------------------------------------------------------------------*/

#ifndef floaterMotionRestraint_H
#define floaterMotionRestraint_H

#include "Time.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "vector.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class floaterMotion;


/*---------------------------------------------------------------------------*\
                Class floaterMotionRestraint Declaration
\*---------------------------------------------------------------------------*/

class floaterMotionRestraint
{

protected:

    // Protected data

        //- Name of the restraint
        word name_;

        //- Restraint model specific coefficient dictionary
        dictionary rBMRCoeffs_;


public:

    //- Runtime type information
    TypeName("floaterMotionRestraint");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            floaterMotionRestraint,
            dictionary,
            (const word& name, const dictionary& rBMRDict),
            (name, rBMRDict)
        );


    // Constructors

        //- Construct from the rBMRDict dictionary and Time
        floaterMotionRestraint
        (
            const word& name,
            const dictionary& rBMRDict
        );

        //- Construct and return a clone
        virtual autoPtr<floaterMotionRestraint> clone() const = 0;


    // Selectors

        //- Select constructed from the rBMRDict dictionary and Time
        static autoPtr<floaterMotionRestraint> New
        (
            const word& name,
            const dictionary& rBMRDict
        );


    //- Destructor
    virtual ~floaterMotionRestraint();


    // Member Functions

        //- Return the name
        const word& name() const
        {
            return name_;
        }

        //- Calculate the restraint position, force and moment.
        //  Global reference frame vectors.
        virtual void restrain
        (
            const floaterMotion& motion,
            vector& restraintPosition,
            vector& restraintForce,
            vector& restraintMoment
        ) const = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& rBMRDict);

        // Access

            // Return access to rBMRCoeffs
            inline const dictionary& coeffDict() const
            {
                return rBMRCoeffs_;
            }

        //- Write
        virtual void write(Ostream&) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
