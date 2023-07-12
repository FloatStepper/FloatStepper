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
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

Class
    Foam::functionObjects::sixDoFRigidBodyStateWriter

Group

Description
    Writes the 6-DoF motion state.

    Example of function object specification:
    \verbatim
    bodyMotion
    {
        type           sixDoFRigidBodyStateWriter;
        libs           (sixDoFRigidBodyMotion);
        angleFormat    degrees;
    }
    \endverbatim

Usage
    \table
        Property     | Description                  | Required | Default value
        type         | type name: sixDoFRigidBodyStateWriter    | yes |
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    sixDoFRigidBodyStateWriter.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyStateWriter_H
#define sixDoFRigidBodyStateWriter_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class sixDoFRigidBodyStateWriter Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyStateWriter
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

        // Old time used to calculate body acceleration (private data so cannot
        // be retrieved).
        scalar t0_;

        // Old body velocity used to calculate acceleration.
        vector v0_;

        // Old body angular velocity used to calculate body angular
        // acceleration.
        vector omega0_;


    // Private Member Functions

        //- No copy construct
        sixDoFRigidBodyStateWriter(const sixDoFRigidBodyStateWriter&) = delete;

        //- No copy assignment
        void operator=(const sixDoFRigidBodyStateWriter&) = delete;


protected:

    // Protected Member Functions

        //- Overloaded writeFileHeader from writeFile
        virtual void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyStateWriter");


    // Constructors

        //- Construct from Time and dictionary
        sixDoFRigidBodyStateWriter
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sixDoFRigidBodyStateWriter() = default;


    // Member Functions

        //- Read the sixDoFRigidBodyStateWriter data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the sixDoFRigidBodyStateWriter
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
