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
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "floaterMotion.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::floaterMotion::read(const dictionary& dict)
{
    dict.readEntry("mass", mass_);
    dict.readEntry("momentOfInertia", momentOfInertia_);

    restraints_.clear();
    addRestraints(dict);

    //constraints_.clear();
    //addConstraints(dict);

    return true;
}


void Foam::floaterMotion::write(Ostream& os) const
{
    motionState_.write(os);

    os.writeEntry("centreOfMass", initialCentreOfMass_);
    os.writeEntry("initialOrientation", initialQ_);
    os.writeEntry("mass", mass_);
    os.writeEntry("momentOfInertia", momentOfInertia_);

    if (!restraints_.empty())
    {
        os.beginBlock("restraints");

        forAll(restraints_, rI)
        {
            const word& restraintType(restraints_[rI].type());

            os.beginBlock(restraints_[rI].name());

            os.writeEntry("floaterMotionRestraint", restraintType);

            restraints_[rI].write(os);

            os.endBlock();
        }

        os.endBlock();
    }
/*
    if (!constraints_.empty())
    {
        os.beginBlock("constraints");

        forAll(constraints_, rI)
        {
            const word& constraintType(constraints_[rI].type());

            os.beginBlock(constraints_[rI].name());

            os.writeEntry("floaterMotionConstraint", constraintType);

            constraints_[rI].floaterMotionConstraint::write(os);

            constraints_[rI].write(os);

            os.endBlock();
        }

        os.endBlock();
    }
*/
}


// ************************************************************************* //
