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

Class
    Foam::slipVelocityFvPatchVectorField

Group
    grpWallBoundaryConditions grpGenericBoundaryConditions

Description
    This boundary condition provides a rigid body motion boundary condition
    with slip where the velocity field is set to rigid body motion with velocity,
    angular velocity and centre of rotation.

    Copyright, Johan Roenby, STROMNING, 2023

SourceFiles
    slipVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef slipVelocityFvPatchVectorField_H
#define slipVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class slipVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class slipVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Current time index (used for updating)
        label curTimeIndex_;

public:

    //- Runtime type information
    TypeName("slipVelocity");


    // Constructors

        //- Construct from patch and internal field
        slipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        slipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given slipVelocityFvPatchVectorField
        //  onto a new patch
        slipVelocityFvPatchVectorField
        (
            const slipVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        slipVelocityFvPatchVectorField
        (
            const slipVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new slipVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        slipVelocityFvPatchVectorField
        (
            const slipVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
             const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new slipVelocityFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        // Access functions

        //- Update the patch vector field value from the given vector field
        virtual void updateCoeffs
        (
            const vector& CoR,
            const vector& v,
            const vector& omega
        );

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
