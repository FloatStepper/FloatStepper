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
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2018 OpenCFD Ltd.
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
    Foam::functionObjects::forceContributions

Group
    grpForcesFunctionObjects

Description
    Derived from forces class with functions to write out pressure, viscous and
    porous force and torque constributions separately.

SourceFiles
    forceContributions.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_forceContributions_H
#define functionObjects_forceContributions_H

#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class forceContributions Declaration
\*---------------------------------------------------------------------------*/

class forceContributions
:
    public forces
{

protected:

public:

    //- Runtime type information
    TypeName("forceContributions");


    // Constructors

        //- Construct from Time and dictionary
        forceContributions
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        forceContributions
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );


    //- Destructor
    virtual ~forceContributions() = default;


    // Member Functions

        //- Return pressure part of force
        virtual vector pressureForce() const;

        //- Return pressure part of moment
        virtual vector pressureMoment() const;

        //- Return viscous part of force
        virtual vector viscousForce() const;

        //- Return viscous part of moment
        virtual vector viscousMoment() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
