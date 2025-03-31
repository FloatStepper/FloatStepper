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
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
    Copyright (C) 2023 Johan Roenby
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
    Foam::rigidBodyVelocityPotentialFvPatchScalarField

Group

Description
    This boundary condition sets the normal gradient of the velocity potential
    according to the specified translational and rotational motion of the body.

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type                rigidBodyVelocityPotential;
        centreOfRotation    (0 0 0);
        linearVelocity      (1 1 0);
        rotationalVelocity  (0 0 1);
    }
    \endverbatim

See also
    Foam::fixedGradientFvPatchField

SourceFiles
    rigidBodyVelocityPotentialFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBodyVelocityPotentialFvPatchScalarField_H
#define rigidBodyVelocityPotentialFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class rigidBodyVelocityPotentialFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class rigidBodyVelocityPotentialFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Current time index (used for updating)
        label curTimeIndex_;

        point CoR_;

        vector v_;

        vector omega_;


public:

    //- Runtime type information
    TypeName("rigidBodyVelocityPotential");


    // Constructors

        //- Construct from patch and internal field
        rigidBodyVelocityPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rigidBodyVelocityPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given rigidBodyVelocityPotentialFvPatchScalarField onto
        //  a new patch
        rigidBodyVelocityPotentialFvPatchScalarField
        (
            const rigidBodyVelocityPotentialFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        rigidBodyVelocityPotentialFvPatchScalarField
        (
            const rigidBodyVelocityPotentialFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new rigidBodyVelocityPotentialFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        rigidBodyVelocityPotentialFvPatchScalarField
        (
            const rigidBodyVelocityPotentialFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new rigidBodyVelocityPotentialFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the patch potential gradient field from the given snGradPhi
        virtual void updateSnGrad(const scalarField& snGradPhi);

        //- Update the patch pressure gradient field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
