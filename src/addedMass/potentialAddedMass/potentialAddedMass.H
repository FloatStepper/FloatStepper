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
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2023 Johan Roenby
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
    potentialAddedMass

Description

SourceFiles
    potentialAddedMass.C

\*---------------------------------------------------------------------------*/

#ifndef potentialAddedMass_H
#define potentialAddedMass_H

#include "addedMass.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class sphereImplicitFunction Declaration
\*---------------------------------------------------------------------------*/

class potentialAddedMass
:
    public addedMass
{

    volScalarField Pot_;

public:

    //- Runtime type information
    TypeName("potentialAddedMass");


    // Constructors

        //- Construct from components
        potentialAddedMass
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~potentialAddedMass() = default;

    virtual SquareMatrix<scalar> computeAddedMass
    (
        FixedList<bool,6> enabledDoF,
        const vector& CoR,
        const tensor& Q,
        const wordRes patches
    );


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
