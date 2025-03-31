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

Class
    Foam::functionObjects::floaterStateWriter

Group

Description
    Writes the 6-DoF motion state for floater motion calculated with
    FloatStepper.

    Example of function object specification:
    \verbatim
    myFloater
    {
        type           floaterStateWriter;
        libs           (floaterStateWriter);
        classicWriteFormat no; //If no, orientation matrix is written out
        angleFormat    degrees;
        writeControl    timeStep;
        writeInterval   1;
        bodyName        circle;
    }
    \endverbatim

Usage
    \table
        Property     | Description                  | Required | Default value
        type         | type name: floaterStateWriter    | yes |
        classicWriteFormat | Use sixDoFRigidBodyStateWrite format | no | yes
        angleFormat | If classicWriteFormat, degrees or radians? | no | radians
        bodyName | Name of body in dynamicMeshDict | no | dynamicMeshDict
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    floaterStateWriter.C

\*---------------------------------------------------------------------------*/

#ifndef floaterStateWriter_H
#define floaterStateWriter_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class floaterStateWriter Declaration
\*---------------------------------------------------------------------------*/

class floaterStateWriter
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private data

        //- The motion solver name
        const word bodyName_;

        enum class angleTypes
        {
            RADIANS,        //!< Radians
            DEGREES         //!< Degrees
        };

        //- Angle type names
        static const Enum<angleTypes> angleTypeNames_;

        //- Angle format
        angleTypes angleFormat_;

        //- Choose classic or raw write format
        Switch classicWriteFormat_;


    // Private Member Functions

        //- No copy construct
        floaterStateWriter(const floaterStateWriter&) = delete;

        //- No copy assignment
        void operator=(const floaterStateWriter&) = delete;


protected:

    // Protected Member Functions

        //- Overloaded writeFileHeader from writeFile
        virtual void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("floaterStateWriter");


    // Constructors

        //- Construct from Time and dictionary
        floaterStateWriter
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~floaterStateWriter() = default;


    // Member Functions

        //- Read the floaterStateWriter data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the floaterStateWriter
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
