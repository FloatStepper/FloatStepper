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
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "floaterMotionState.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::floaterMotionState::write(dictionary& dict) const
{
    dict.add("centreOfRotation", centreOfRotation_);
    dict.add("orientation", Q_);
    dict.add("velocity", v_);
    dict.add("acceleration", a_);
    dict.add("angularVelocity", omega_);
    dict.add("angularAcceleration", domegadt_);
}


void Foam::floaterMotionState::write(Ostream& os) const
{
    os.writeEntry("centreOfRotation", centreOfRotation_);
    os.writeEntry("orientation", Q_);
    os.writeEntry("velocity", v_);
    os.writeEntry("acceleration", a_);
    os.writeEntry("angularVelocity", omega_);
    os.writeEntry("angularAcceleration", domegadt_);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is, floaterMotionState& rBMS
)
{
    is  >> rBMS.centreOfRotation_
        >> rBMS.Q_
        >> rBMS.v_
        >> rBMS.a_
        >> rBMS.omega_
        >> rBMS.domegadt_;

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const floaterMotionState& rBMS
)
{
    os  << token::SPACE << rBMS.centreOfRotation()
        << token::SPACE << rBMS.Q()
        << token::SPACE << rBMS.v()
        << token::SPACE << rBMS.a()
        << token::SPACE << rBMS.omega()
        << token::SPACE << rBMS.domegadt();

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
