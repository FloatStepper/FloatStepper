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
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "oscillatingForce.H"
#include "addToRunTimeSelectionTable.H"
#include "floaterMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace floaterMotionRestraints
{
    defineTypeNameAndDebug(oscillatingForce, 0);

    addToRunTimeSelectionTable
    (
        floaterMotionRestraint,
        oscillatingForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::floaterMotionRestraints::oscillatingForce::oscillatingForce
(
    const word& name,
    const dictionary& rBMRDict
)
:
    floaterMotionRestraint(name, rBMRDict),
    direction_(),
    period_(),
    amplitude_()
{
    read(rBMRDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::floaterMotionRestraints::oscillatingForce::~oscillatingForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::floaterMotionRestraints::oscillatingForce::restrain
(
    const floaterMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    const scalar time = motion.time().value();

    const scalar omega = Foam::constant::mathematical::twoPi/period_;

    restraintPosition = motion.centreOfRotation();

    restraintForce = amplitude_*(direction_/mag(direction_))
        *Foam::sin(omega*(time - time0_/period_));

    restraintMoment = Zero;

    Info<< "Oscillating force: " << restraintForce << endl;
}


bool Foam::floaterMotionRestraints::oscillatingForce::read
(
    const dictionary& rBMRDict
)
{
    floaterMotionRestraint::read(rBMRDict);

    rBMRCoeffs_.readEntry("direction", direction_);
    rBMRCoeffs_.readEntry("period", period_);
    rBMRCoeffs_.readEntry("amplitude", amplitude_);
    rBMRCoeffs_.readEntry("zeroTime", time0_);

    return true;
}


void Foam::floaterMotionRestraints::oscillatingForce::write
(
    Ostream& os
) const
{
    os.writeEntry("direction", direction_);
    os.writeEntry("period", period_);
    os.writeEntry("amplitude", amplitude_);
    os.writeEntry("zeroTime", time0_);
}

// ************************************************************************* //
