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
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-----------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2023 Johan Roenby
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
    pisoAddedMass

Description

SourceFiles
    pisoAddedMass.C

\*---------------------------------------------------------------------------*/

#ifndef pisoAddedMass_H
#define pisoAddedMass_H

#include "addedMass.H"
#include "forceContributions.H"
#include "slipVelocityFvPatchVectorField.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class sphereImplicitFunction Declaration
\*---------------------------------------------------------------------------*/

class pisoAddedMass
:
    public addedMass
{

//    volScalarField Pot_;
    const fvMesh& mesh_;

public:

    //- Runtime type information
    TypeName("pisoAddedMass");


    // Constructors

        //- Construct from components
        pisoAddedMass
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~pisoAddedMass() = default;

    virtual SquareMatrix<scalar> computeAddedMass
    (
        FixedList<bool,6> enabledDoF,
        const vector& CoR,
        const tensor& Q,
        const wordRes patches
    );


    void calcPressureForceAndTorque
    (
        const fvMesh& mesh,
        const vector& v,
        const vector& omega,
        const point& CoR,
        const wordRes& patches,
        vector& force,
        vector& torque
    );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
