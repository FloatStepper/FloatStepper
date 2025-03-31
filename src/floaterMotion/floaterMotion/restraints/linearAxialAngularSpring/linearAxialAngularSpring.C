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

#include "linearAxialAngularSpring.H"
#include "addToRunTimeSelectionTable.H"
#include "floaterMotion.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace floaterMotionRestraints
{
    defineTypeNameAndDebug(linearAxialAngularSpring, 0);

    addToRunTimeSelectionTable
    (
        floaterMotionRestraint,
        linearAxialAngularSpring,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::floaterMotionRestraints::linearAxialAngularSpring::
linearAxialAngularSpring
(
    const word& name,
    const dictionary& rBMRDict
)
:
    floaterMotionRestraint(name, rBMRDict),
    refQ_(),
    axis_(),
    stiffness_(),
    damping_()
{
    read(rBMRDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::floaterMotionRestraints::linearAxialAngularSpring::
~linearAxialAngularSpring()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void
Foam::floaterMotionRestraints::linearAxialAngularSpring::restrain
(
    const floaterMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    vector refDir = rotationTensor(vector(1, 0 ,0), axis_) & vector(0, 1, 0);

    vector oldDir = refQ_ & refDir;
    vector newDir = motion.orientation() & refDir;

    if (mag(oldDir & axis_) > 0.95 || mag(newDir & axis_) > 0.95)
    {
        // Directions getting close to the axis, change reference

        refDir = rotationTensor(vector(1, 0 ,0), axis_) & vector(0, 0, 1);
        oldDir = refQ_ & refDir;
        newDir = motion.orientation() & refDir;
    }

    // Removing any axis component from oldDir and newDir and normalising
    oldDir -= (axis_ & oldDir)*axis_;
    oldDir /= (mag(oldDir) + VSMALL);

    newDir -= (axis_ & newDir)*axis_;
    newDir /= (mag(newDir) + VSMALL);

    scalar theta = mag(acos(min(oldDir & newDir, 1.0)));

    // Temporary axis with sign information.
    vector a = (oldDir ^ newDir);

    // Remove any component that is not along axis that may creep in
    a = (a & axis_)*axis_;

    scalar magA = mag(a);

    if (magA > VSMALL)
    {
        a /= magA;
    }
    else
    {
        a = Zero;
    }

    // Damping of along axis angular velocity only
    restraintMoment = -stiffness_*theta*a - damping_*(motion.omega() & a)*a;

    restraintForce = Zero;

    // Not needed to be altered as restraintForce is zero, but set to
    // centreOfRotation to be sure of no spurious moment
    restraintPosition = motion.centreOfRotation();

    Info<< " angle " << theta*sign(a & axis_) << " moment "
        << restraintMoment << endl;
}


bool Foam::floaterMotionRestraints::linearAxialAngularSpring::read
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

    rBMRCoeffs_.readEntry("axis", axis_);

    const scalar magAxis(mag(axis_));

    if (magAxis > VSMALL)
    {
        axis_ /= magAxis;
    }
    else
    {
        FatalErrorInFunction
            << "axis has zero length"
            << abort(FatalError);
    }

    rBMRCoeffs_.readEntry("stiffness", stiffness_);
    rBMRCoeffs_.readEntry("damping", damping_);

    return true;
}


void Foam::floaterMotionRestraints::linearAxialAngularSpring::write
(
    Ostream& os
) const
{
    os.writeEntry("referenceOrientation", refQ_);
    os.writeEntry("axis", axis_);
    os.writeEntry("stiffness", stiffness_);
    os.writeEntry("damping", damping_);
}

// ************************************************************************* //
