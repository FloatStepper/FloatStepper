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

\*---------------------------------------------------------------------------*/

#include "linearSpring.H"
#include "addToRunTimeSelectionTable.H"
#include "floaterMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace floaterMotionRestraints
{
    defineTypeNameAndDebug(linearSpring, 0);

    addToRunTimeSelectionTable
    (
        floaterMotionRestraint,
        linearSpring,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::floaterMotionRestraints::linearSpring::linearSpring
(
    const word& name,
    const dictionary& rBMRDict
)
:
    floaterMotionRestraint(name, rBMRDict),
    anchor_(),
    refAttachmentPt_(),
    stiffness_(),
    damping_(),
    restLength_()
{
    read(rBMRDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::floaterMotionRestraints::linearSpring::~linearSpring()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::floaterMotionRestraints::linearSpring::restrain
(
    const floaterMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    restraintPosition = motion.transform(refAttachmentPt_);

    vector r = restraintPosition - anchor_;

    scalar magR = mag(r);
    r /= (magR + VSMALL);

    vector v = motion.velocity(restraintPosition);

//    restraintForce = -stiffness_*(magR - restLength_)*r - damping_*(r & v)*r;
    restraintForce = -stiffness_*(magR - restLength_)*r - damping_*v;

    restraintMoment = Zero;

    Info<< " attachmentPt - anchor " << r*magR
        << " spring length " << magR
        << " force " << restraintForce
        << endl;
}


bool Foam::floaterMotionRestraints::linearSpring::read
(
    const dictionary& rBMRDict
)
{
    floaterMotionRestraint::read(rBMRDict);

    rBMRCoeffs_.readEntry("anchor", anchor_);
    rBMRCoeffs_.readEntry("refAttachmentPt", refAttachmentPt_);
    rBMRCoeffs_.readEntry("stiffness", stiffness_);
    rBMRCoeffs_.readEntry("damping", damping_);
    rBMRCoeffs_.readEntry("restLength", restLength_);

    return true;
}


void Foam::floaterMotionRestraints::linearSpring::write
(
    Ostream& os
) const
{
    os.writeEntry("anchor", anchor_);
    os.writeEntry("refAttachmentPt", refAttachmentPt_);
    os.writeEntry("stiffness", stiffness_);
    os.writeEntry("damping", damping_);
    os.writeEntry("restLength", restLength_);
}

// ************************************************************************* //
