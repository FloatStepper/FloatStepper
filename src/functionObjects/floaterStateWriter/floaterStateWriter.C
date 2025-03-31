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
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "floaterStateWriter.H"
#include "floaterMotionSolver.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(floaterStateWriter, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        floaterStateWriter,
        dictionary
    );
}
}


const Foam::Enum
<
    Foam::functionObjects::floaterStateWriter::angleTypes
>
Foam::functionObjects::floaterStateWriter::angleTypeNames_
({
    { angleTypes::RADIANS, "radians" },
    { angleTypes::DEGREES, "degrees" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::floaterStateWriter::writeFileHeader(Ostream& os)
{
    if (classicWriteFormat_)
    {
        writeHeader(os, "Motion State");
        writeHeaderValue(os, "Angle Units", angleTypeNames_[angleFormat_]);
        writeCommented(os, "Time");

        os  << tab
            << "centreOfRotation" << tab
            << "centreOfMass" << tab
            << "rotation" << tab
            << "velocity" << tab
            << "angularVelocity" << endl;
    }
    else
    {
        writeHeader(os, "Rigid body motion state for body " + bodyName_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::floaterStateWriter::floaterStateWriter
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    bodyName_(dict.getOrDefault<word>("bodyName", "dynamicMeshDict")),
    angleFormat_(angleTypes::RADIANS),
    classicWriteFormat_(dict.getOrDefault<Switch>("classicWriteFormat", true))
{
    read(dict);
    writeFileHeader(file());
    write();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::floaterStateWriter::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {

        angleFormat_ =
            angleTypeNames_.lookupOrDefault
            (
                "angleFormat",
                dict,
                angleTypes::RADIANS
            );

        return true;
    }

    return false;
}


bool Foam::functionObjects::floaterStateWriter::execute()
{
    return true;
}


bool Foam::functionObjects::floaterStateWriter::write()
{

    const floaterMotionSolver& bodySolver =
        mesh_.lookupObject<floaterMotionSolver>(bodyName_);
    const floaterMotion& motion(bodySolver.motion());
    const floaterMotionState& state(motion.state());

    if (classicWriteFormat_)
    {

        vector rotationAngle
        (
            quaternion(motion.orientation()).eulerAngles(quaternion::XYZ)
        );

        vector angularVelocity(motion.omega());

        switch (angleFormat_)
        {
            case angleTypes::RADIANS:
            {
                // Nothing to do - already in radians
                break;
            }
            case angleTypes::DEGREES:
            {
                rotationAngle.x() = radToDeg(rotationAngle.x());
                rotationAngle.y() = radToDeg(rotationAngle.y());
                rotationAngle.z() = radToDeg(rotationAngle.z());

                angularVelocity.x() = radToDeg(angularVelocity.x());
                angularVelocity.y() = radToDeg(angularVelocity.y());
                angularVelocity.z() = radToDeg(angularVelocity.z());
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unhandled enumeration " << angleTypeNames_[angleFormat_]
                    << abort(FatalError);
            }
        }

        writeCurrentTime(file());
        file()
            << tab
            << motion.centreOfRotation()  << tab
            << motion.centreOfMass()  << tab
            << rotationAngle  << tab
            << motion.v()  << tab
            << angularVelocity << endl;
    }
    else
    {
        file() << "Time = ";
        writeCurrentTime(file());
        file()
            << endl
            << "Centre of rotation  : " << state.centreOfRotation() << endl
            << "Velocity            : " << state.v() << endl
            << "Acceleration        : " << state.a() << endl
            << "Orientation         : " << state.Q() << endl
            << "Angular velocity    : " << state.omega() << endl
            << "Angular acceleration: " << state.domegadt() << endl
            << endl;
    }

    return true;
}


// ************************************************************************* //
