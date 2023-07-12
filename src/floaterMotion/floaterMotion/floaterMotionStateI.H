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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::point&
Foam::floaterMotionState::centreOfRotation() const
{
    return centreOfRotation_;
}


inline const Foam::tensor& Foam::floaterMotionState::Q() const
{
    return Q_;
}


inline const Foam::vector& Foam::floaterMotionState::v() const
{
    return v_;
}


inline const Foam::vector& Foam::floaterMotionState::a() const
{
    return a_;
}


inline const Foam::vector& Foam::floaterMotionState::omega() const
{
    return omega_;
}


inline const Foam::vector& Foam::floaterMotionState::domegadt() const
{
    return domegadt_;
}


inline Foam::point& Foam::floaterMotionState::centreOfRotation()
{
    return centreOfRotation_;
}


inline Foam::tensor& Foam::floaterMotionState::Q()
{
    return Q_;
}


inline Foam::vector& Foam::floaterMotionState::v()
{
    return v_;
}


inline Foam::vector& Foam::floaterMotionState::a()
{
    return a_;
}


inline Foam::vector& Foam::floaterMotionState::omega()
{
    return omega_;
}


inline Foam::vector& Foam::floaterMotionState::domegadt()
{
    return domegadt_;
}


// ************************************************************************* //
