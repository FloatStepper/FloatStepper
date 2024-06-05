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
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "floaterMotionState.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::floaterMotionState::floaterMotionState()
:
    centreOfRotation_(Zero),
    Q_(I),
    v_(Zero),
    a_(Zero),
    omega_(Zero),
    domegadt_(Zero)
{}

/*
Foam::floaterMotionState::floaterMotionState
(
    const dictionary& dict
)
:
    centreOfRotation_
    (
        dict.getOrDefault
        (
            "centreOfRotation",
            dict.getOrDefault("centreOfMass", vector::zero)
        )
    ),
    Q_(dict.getOrDefault("orientation", tensor::I)),
    v_(dict.getOrDefault("velocity", vector::zero)),
    a_(dict.getOrDefault("acceleration", vector::zero)),
    omega_(dict.getOrDefault("omega", vector::zero)),
    domegadt_(dict.getOrDefault("domegadt", vector::zero))
{}
*/

// New reader enforcing explicit specification of all motion parameters
Foam::floaterMotionState::floaterMotionState
(
    const dictionary& dict
)
:
    centreOfRotation_(dict.get<point>("centreOfRotation")),
    Q_(dict.get<tensor>("orientation")),
    v_(dict.get<vector>("velocity")),
    a_(dict.get<vector>("acceleration")),
    omega_(dict.get<vector>("omega")),
    domegadt_(dict.get<vector>("domegadt"))
{}


Foam::floaterMotionState::floaterMotionState
(
    const floaterMotionState& rBMS
)
:
    centreOfRotation_(rBMS.centreOfRotation()),
    Q_(rBMS.Q()),
    v_(rBMS.v()),
    a_(rBMS.a()),
    omega_(rBMS.omega()),
    domegadt_(rBMS.domegadt())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::floaterMotionState::~floaterMotionState()
{}


// ************************************************************************* //
