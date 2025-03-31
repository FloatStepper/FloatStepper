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
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


inline const Foam::PtrList<Foam::floaterMotionRestraint>&
Foam::floaterMotion::restraints() const
{
    return restraints_;
}

/*
inline const Foam::PtrList<Foam::floaterMotionConstraint>&
Foam::floaterMotion::constraints() const
{
    return constraints_;
}
*/

inline const Foam::point&
Foam::floaterMotion::initialCentreOfRotation() const
{
    return initialCentreOfRotation_;
}


inline const Foam::tensor&
Foam::floaterMotion::initialQ() const
{
    return initialQ_;
}


inline const Foam::tensor& Foam::floaterMotion::Q() const
{
    return motionState_.Q();
}


inline const Foam::vector& Foam::floaterMotion::v() const
{
    return motionState_.v();
}


inline Foam::point& Foam::floaterMotion::initialCentreOfRotation()
{
    return initialCentreOfRotation_;
}


inline Foam::tensor& Foam::floaterMotion::initialQ()
{
    return initialQ_;
}


inline Foam::tensor& Foam::floaterMotion::Q()
{
    return motionState_.Q();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


inline Foam::word Foam::floaterMotion::name() const
{
    return name_;
}


inline Foam::scalar Foam::floaterMotion::mass() const
{
    return mass_;
}


inline const Foam::symmTensor&
Foam::floaterMotion::momentOfInertia() const
{
    return momentOfInertia_;
}


inline const Foam::floaterMotionState&
Foam::floaterMotion::state() const
{
    return motionState_;
}


inline const Foam::point& Foam::floaterMotion::centreOfRotation() const
{
    return motionState_.centreOfRotation();
}


inline const Foam::point&
Foam::floaterMotion::initialCentreOfMass() const
{
    return initialCentreOfMass_;
}


inline Foam::point Foam::floaterMotion::centreOfMass() const
{
    return transform(initialCentreOfMass_);
}


inline Foam::vector Foam::floaterMotion::momentArm() const
{
    return centreOfMass() - motionState_.centreOfRotation();
}


inline const Foam::Time& Foam::floaterMotion::time() const
{
    return time_;
}


inline const Foam::vector&
Foam::floaterMotion::restraintForce() const
{
    return restraintForce_;
}


inline const Foam::vector&
Foam::floaterMotion::restraintTorque() const
{
    return restraintTorque_;
}


inline const Foam::scalarSquareMatrix&
Foam::floaterMotion::Madd() const
{
    return Madd_;
}


inline const Foam::vector&
Foam::floaterMotion::F0() const
{
    return F0_;
}


inline const Foam::vector&
Foam::floaterMotion::tau0() const
{
    return tau0_;
}


inline const Foam::tensor&
Foam::floaterMotion::orientation() const
{
    return Q();
}


inline const Foam::vector& Foam::floaterMotion::omega() const
{
    return  motionState_.omega();
}


inline void Foam::floaterMotion::storeOldState()
{
    motionState0_ = motionState_;
}


inline Foam::point& Foam::floaterMotion::centreOfRotation()
{
    return motionState_.centreOfRotation();
}


inline void Foam::floaterMotion::setVelocity
(
    const vector& V
)
{
    motionState_.v() = V;
}

inline void Foam::floaterMotion::setAngularVelocity
(
    const vector& omega
)
{
    motionState_.omega() = omega;
}


inline void Foam::floaterMotion::setCentreOfRotation
(
    const vector& CoR
)
{
    motionState_.centreOfRotation() = CoR;
}


inline void Foam::floaterMotion::setOrientation
(
    const tensor& Q
)
{
    motionState_.Q() = Q;
}


inline void Foam::floaterMotion::setAcceleration
(
    const vector& a
)
{
    motionState_.a() = a;
}


inline void Foam::floaterMotion::setAngularAcceleration
(
    const vector& domegadt
)
{
    motionState_.domegadt() = domegadt;
}


inline void Foam::floaterMotion::setState
(
    const floaterMotionState& state
)
{
    motionState_ = state;
}


inline Foam::point Foam::floaterMotion::velocity
(
    const point& pt
) const
{
    return (omega() ^ (pt - centreOfRotation())) + v();
}


inline Foam::point Foam::floaterMotion::transform
(
    const point& initialPoint
) const
{
    return
    (
        centreOfRotation()
      + (Q() & initialQ_.T() & (initialPoint - initialCentreOfRotation_))
    );
}


// ************************************************************************* //
