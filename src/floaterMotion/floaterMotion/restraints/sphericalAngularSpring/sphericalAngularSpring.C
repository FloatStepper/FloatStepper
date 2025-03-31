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
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "sphericalAngularSpring.H"
#include "addToRunTimeSelectionTable.H"
#include "floaterMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace floaterMotionRestraints
{
    defineTypeNameAndDebug(sphericalAngularSpring, 0);

    addToRunTimeSelectionTable
    (
        floaterMotionRestraint,
        sphericalAngularSpring,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::floaterMotionRestraints::sphericalAngularSpring::
sphericalAngularSpring
(
    const word& name,
    const dictionary& rBMRDict
)
:
    floaterMotionRestraint(name, rBMRDict),
    refQ_(),
    stiffness_(),
    damping_()
{
    read(rBMRDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::floaterMotionRestraints::sphericalAngularSpring::
~sphericalAngularSpring()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::floaterMotionRestraints::sphericalAngularSpring::restrain
(
    const floaterMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    restraintMoment = Zero;

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        vector axis = Zero;
        axis[cmpt] = 1;

        vector refDir = Zero;
        refDir[(cmpt + 1) % 3] = 1;

        vector newDir = motion.orientation() & refDir;

        axis = (refQ_ & axis);
        refDir = (refQ_ & refDir);
        newDir -= (axis & newDir)*axis;

        restraintMoment += -stiffness_*(refDir ^ newDir);
    }

    restraintMoment += -damping_*motion.omega();

    restraintForce = Zero;

    // Not needed to be altered as restraintForce is zero, but set to
    // centreOfRotation to be sure of no spurious moment
    restraintPosition = motion.centreOfRotation();

    Info<< " moment " << restraintMoment << endl;
}


bool Foam::floaterMotionRestraints::sphericalAngularSpring::read
(
    const dictionary& rBMRDict
)
{
    floaterMotionRestraint::read(rBMRDict);

    refQ_ = rBMRCoeffs_.getOrDefault<tensor>("referenceOrientation", I);

    if (mag(mag(refQ_) - sqrt(3.0)) > ROOTSMALL)
    {
        FatalErrorInFunction
            << "referenceOrientation " << refQ_ << " is not a rotation tensor. "
            << "mag(referenceOrientation) - sqrt(3) = "
            << mag(refQ_) - sqrt(3.0) << nl
            << exit(FatalError);
    }

    rBMRCoeffs_.readEntry("stiffness", stiffness_);
    rBMRCoeffs_.readEntry("damping", damping_);

    return true;
}


void Foam::floaterMotionRestraints::sphericalAngularSpring::write
(
    Ostream& os
) const
{
    os.writeEntry("referenceOrientation", refQ_);
    os.writeEntry("stiffness", stiffness_);
    os.writeEntry("damping", damping_);
}


// ************************************************************************* //
