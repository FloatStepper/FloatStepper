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
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "floaterMotion.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::floaterMotion::read(const dictionary& dict)
{
    dict.readEntry("mass", mass_);
    dict.readEntry("momentOfInertia", momentOfInertia_);

    Info << "Warning: momentOfInertia read by floaterMotion::read() but " 
        << "momentOfInertia not transformed!" << endl;
    Info << "Please make sure it also reads momentOfInertiaRefPoint and "
        << "momentOfInertiaAxes and transform accordingly." << endl;

    restraints_.clear();
    addRestraints(dict);

    //constraints_.clear();
    //addConstraints(dict);

    return true;
}


void Foam::floaterMotion::write(dictionary& dict) const
{
    motionState_.write(dict);

    dict.add("mass", mass_);
    dict.add("centreOfMass", centreOfMass());
    dict.add("momentOfInertia", momentOfInertia_);
//    dict.add("momentOfInertiaRefPoint", momentOfInertiaRefPoint_);
//    dict.add("momentOfInertiaAxes", momentOfInertiaAxes_);
    dict.add("initialCentreOfRotation", initialCentreOfRotation_);
    dict.add("initialOrientation", initialQ_);
    dict.add("initialCentreOfMass", initialCentreOfMass_);
    dict.add("addedMass", Madd_);
    dict.add("F0", F0_);
    dict.add("tau0", tau0_);

}


void Foam::floaterMotion::write(Ostream& os) const
{
    motionState_.write(os);

    os.writeEntry("mass", mass_);
    os.writeEntry("centreOfMass", centreOfMass());
    os.writeEntry("momentOfInertia", momentOfInertia_);
//    os.writeEntry("momentOfInertiaRefPoint", momentOfInertiaRefPoint_);
//    os.writeEntry("momentOfInertiaAxes", momentOfInertiaAxes_);
    os.writeEntry("initialCentreOfRotation", initialCentreOfRotation_);
    os.writeEntry("initialOrientation", initialQ_);
    os.writeEntry("initialCentreOfMass", initialCentreOfMass_);
    os.writeEntry("addedMass", Madd_);
    os.writeEntry("F0", F0_);
    os.writeEntry("tau0", tau0_);

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
