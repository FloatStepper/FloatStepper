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
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "tensor.H"
#include "cubicEqn.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::tensor::vsType::typeName = "tensor";

template<>
const char* const Foam::tensor::vsType::componentNames[] =
{
    "xx", "xy", "xz",
    "yx", "yy", "yz",
    "zx", "zy", "zz"
};

template<>
const Foam::tensor Foam::tensor::vsType::zero(tensor::uniform(0));

template<>
const Foam::tensor Foam::tensor::vsType::one(tensor::uniform(1));

template<>
const Foam::tensor Foam::tensor::vsType::max(tensor::uniform(VGREAT));

template<>
const Foam::tensor Foam::tensor::vsType::min(tensor::uniform(-VGREAT));

template<>
const Foam::tensor Foam::tensor::vsType::rootMax(tensor::uniform(ROOTVGREAT));

template<>
const Foam::tensor Foam::tensor::vsType::rootMin(tensor::uniform(-ROOTVGREAT));

template<>
const Foam::tensor Foam::tensor::I
(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::vector Foam::eigenValues(const tensor& t)
{
    // Coefficients of the characteristic cubic polynomial (a = 1)
    const scalar b =
      - t.xx() - t.yy() - t.zz();
    const scalar c =
        t.xx()*t.yy() + t.xx()*t.zz() + t.yy()*t.zz()
      - t.xy()*t.yx() - t.yz()*t.zy() - t.zx()*t.xz();
    const scalar d =
      - t.xx()*t.yy()*t.zz()
      - t.xy()*t.yz()*t.zx() - t.xz()*t.zy()*t.yx()
      + t.xx()*t.yz()*t.zy() + t.yy()*t.zx()*t.xz() + t.zz()*t.xy()*t.yx();

    // Solve
    Roots<3> roots = cubicEqn(1, b, c, d).roots();

    // Check the root types
    vector lambda = vector::zero;
    forAll(roots, i)
    {
        switch (roots.type(i))
        {
            case roots::real:
                lambda[i] = roots[i];
                break;
            case roots::complex:
                WarningInFunction
                    << "Complex eigenvalues detected for tensor: " << t
                    << endl;
                lambda[i] = 0;
                break;
            case roots::posInf:
                lambda[i] = VGREAT;
                break;
            case roots::negInf:
                lambda[i] = - VGREAT;
                break;
            case roots::nan:
                FatalErrorInFunction
                    << "Eigenvalue calculation failed for tensor: " << t
                    << exit(FatalError);
        }
    }

    // Sort the eigenvalues into ascending order
    if (lambda.x() > lambda.y())
    {
        Swap(lambda.x(), lambda.y());
    }
    if (lambda.y() > lambda.z())
    {
        Swap(lambda.y(), lambda.z());
    }
    if (lambda.x() > lambda.y())
    {
        Swap(lambda.x(), lambda.y());
    }

    return lambda;
}


Foam::vector Foam::eigenVector
(
    const tensor& T,
    const scalar lambda,
    const vector& direction1,
    const vector& direction2
)
{
    // Construct the linear system for this eigenvalue
    tensor A(T - lambda*I);

    // Determinants of the 2x2 sub-matrices used to find the eigenvectors
    scalar sd0, sd1, sd2;
    scalar magSd0, magSd1, magSd2;

    // Sub-determinants for a unique eivenvalue
    sd0 = A.yy()*A.zz() - A.yz()*A.zy();
    sd1 = A.zz()*A.xx() - A.zx()*A.xz();
    sd2 = A.xx()*A.yy() - A.xy()*A.yx();
    magSd0 = mag(sd0);
    magSd1 = mag(sd1);
    magSd2 = mag(sd2);

    // Evaluate the eigenvector using the largest sub-determinant
    if (magSd0 >= magSd1 && magSd0 >= magSd2 && magSd0 > SMALL)
    {
        vector ev
        (
            1,
            (A.yz()*A.zx() - A.zz()*A.yx())/sd0,
            (A.zy()*A.yx() - A.yy()*A.zx())/sd0
        );

        return ev/mag(ev);
    }
    else if (magSd1 >= magSd2 && magSd1 > SMALL)
    {
        vector ev
        (
            (A.xz()*A.zy() - A.zz()*A.xy())/sd1,
            1,
            (A.zx()*A.xy() - A.xx()*A.zy())/sd1
        );

        return ev/mag(ev);
    }
    else if (magSd2 > SMALL)
    {
        vector ev
        (
            (A.xy()*A.yz() - A.yy()*A.xz())/sd2,
            (A.yx()*A.xz() - A.xx()*A.yz())/sd2,
            1
        );

        return ev/mag(ev);
    }

    // Sub-determinants for a repeated eigenvalue
    sd0 = A.yy()*direction1.z() - A.yz()*direction1.y();
    sd1 = A.zz()*direction1.x() - A.zx()*direction1.z();
    sd2 = A.xx()*direction1.y() - A.xy()*direction1.x();
    magSd0 = mag(sd0);
    magSd1 = mag(sd1);
    magSd2 = mag(sd2);

    // Evaluate the eigenvector using the largest sub-determinant
    if (magSd0 >= magSd1 && magSd0 >= magSd2 && magSd0 > SMALL)
    {
        vector ev
        (
            1,
            (A.yz()*direction1.x() - direction1.z()*A.yx())/sd0,
            (direction1.y()*A.yx() - A.yy()*direction1.x())/sd0
        );

        return ev/mag(ev);
    }
    else if (magSd1 >= magSd2 && magSd1 > SMALL)
    {
        vector ev
        (
            (direction1.z()*A.zy() - A.zz()*direction1.y())/sd1,
            1,
            (A.zx()*direction1.y() - direction1.x()*A.zy())/sd1
        );

        return ev/mag(ev);
    }
    else if (magSd2 > SMALL)
    {
        vector ev
        (
            (A.xy()*direction1.z() - direction1.y()*A.xz())/sd2,
            (direction1.x()*A.xz() - A.xx()*direction1.z())/sd2,
            1
        );

        return ev/mag(ev);
    }

    // Triple eigenvalue
    return direction1^direction2;
}


Foam::tensor Foam::eigenVectors(const tensor& T, const vector& lambdas)
{
    vector Ux(1, 0, 0), Uy(0, 1, 0), Uz(0, 0, 1);

    Ux = eigenVector(T, lambdas.x(), Uy, Uz);
    Uy = eigenVector(T, lambdas.y(), Uz, Ux);
    Uz = eigenVector(T, lambdas.z(), Ux, Uy);

    return tensor(Ux, Uy, Uz);
}


Foam::tensor Foam::eigenVectors(const tensor& T)
{
    const vector lambdas(eigenValues(T));

    return eigenVectors(T, lambdas);
}


Foam::vector Foam::eigenValues(const symmTensor& T)
{
    // Reference:
    //    \verbatim
    //        Smith, Oliver K. (April 1961)
    //           "Eigenvalues of a symmetric 3 × 3 matrix."
    //           Communications of the ACM, 4 (4): 168
    //           doi:10.1145/355578.366316
    //      \endverbatim

    vector eig;

    scalar p1 = sqr(T.xy()) + sqr(T.xz()) + sqr(T.yz());
    if (p1 < SMALL)
    {
        // T is diagonal.
        eig[0] = T.xx();
        eig[1] = T.yy();
        eig[2] = T.zz();
    }
    else
    {
        scalar q = tr(T)/3;
        scalar p2 = sqr(T.xx() - q) + sqr(T.yy() - q) + sqr(T.zz() - q) + 2*p1;
        scalar p = sqrt(p2/6);
        symmTensor B = (1/p)*(T - q*symmTensor::I);
        scalar r = det(B)/2;

        // In exact arithmetic for a symmetric matrix  -1 <= r <= 1
        // but computation error can leave it slightly outside this range.
        scalar phi(0);
        if (r <= -1)
        {
           phi = pi/3;
        }
        else if (r >= 1)
        {
           phi = 0;
        }
        else
        {
           phi = acos(r)/3;
        }

        // the eigenvalues satisfy eig3 <= eig2 <= eig1
        eig[0] = q + 2*p*cos(phi);
        eig[2] = q + 2*p*cos(phi + (2*pi/3));
        eig[1] = 3*q - eig[0] - eig[2];
    }

    return eig;
}


Foam::vector Foam::eigenVector
(
    const symmTensor& T,
    const scalar lambda,
    const vector& direction1,
    const vector& direction2
)
{
    return eigenVector(tensor(T), lambda, direction1, direction2);
}


Foam::tensor Foam::eigenVectors(const symmTensor& T, const vector& lambdas)
{
    return eigenVectors(tensor(T), lambdas);
}


Foam::tensor Foam::eigenVectors(const symmTensor& T)
{
    return eigenVectors(tensor(T));
}


// ************************************************************************* //
