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
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-----------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "addedMassWriter.H"
#include "floaterMeshMotionSolver.H"
//#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"
#include "MatrixTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(addedMassWriter, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        addedMassWriter,
        dictionary
    );
}
}

/*
const Foam::Enum
<
    Foam::functionObjects::addedMassWriter::angleTypes
>
Foam::functionObjects::addedMassWriter::angleTypeNames_
({
    { angleTypes::RADIANS, "radians" },
    { angleTypes::DEGREES, "degrees" },
});
*/

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::addedMassWriter::writeFileHeader(Ostream& os)
{
    writeHeader(os, "Added mass matrix, F0 and tau0");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::addedMassWriter::addedMassWriter
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
//    angleFormat_(angleTypes::RADIANS),
    bodyName_(dict.getOrDefault<word>("bodyName", "dynamicMeshDict")),
    MaddInBodyFrame_(true)
{
    read(dict);
    writeFileHeader(file());
    write();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::addedMassWriter::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        MaddInBodyFrame_ = dict.lookupOrDefault("MaddInBodyFrame", true);

        return true;
    }

    return false;
}


bool Foam::functionObjects::addedMassWriter::execute()
{
    return true;
}


bool Foam::functionObjects::addedMassWriter::write()
{

    const floaterMeshMotionSolver& bodySolver =
        mesh_.lookupObject<floaterMeshMotionSolver>(bodyName_);

    const floaterMotion& motion(bodySolver.motion());
    scalarSquareMatrix Madd = motion.Madd();
    if (!MaddInBodyFrame_)
    {
        const tensor Q = motion.orientation();
        Madd = motion.changeFrame(Madd, Q.T());
    }

    file() << "Time = ";
    writeCurrentTime(file());
    MatrixTools::printMatrix(file(), Madd);
    file()
        << endl
        << "F0               : " << motion.F0() << endl
        << "tau0             : " << motion.tau0() << endl
        << endl;

    return true;
}


// ************************************************************************* //
