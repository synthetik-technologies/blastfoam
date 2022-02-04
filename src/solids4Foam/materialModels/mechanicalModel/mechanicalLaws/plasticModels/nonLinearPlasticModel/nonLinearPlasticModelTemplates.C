/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "nonLinearPlasticModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<template<class> class PatchField, class Mesh>
Foam::tmp<Foam::GeometricField<Foam::scalar, PatchField, Mesh>>
Foam::nonLinearPlasticModel::Ibar
(
    const GeometricField<symmTensor, PatchField, Mesh>& devBEbar
)
{
    // From Simo & Hughes 1998:
    // but this incorrectly results in det(bEbar) =! 1
    //bEbar = (s/mu) + Ibar*I;

    // A method of calculating Ibar to enforce det(bEbar) == 1 is proposed
    // by solving a cubic equation.
    // Rubin and Attia, CALCULATION OF HYPERELASTIC RESPONSE OF FINITELY
    // DEFORMED ELASTIC-VISCOPLASTIC MATERIALS, INTERNATIONAL JOURNAL FOR
    // NUMERICAL METHODS IN ENGINEERING, VOL. 39,309-320(1996)
    // and
    // M. Hollenstein M. Jabareen M. B. Rubin, Modeling a smooth elastic-
    // inelastic transition with a strongly objective numerical integrator
    // needing no iteration, Comput Mech (2013) 52:649–667
    // DOI 10.1007/s00466-013-0838-7

    // Note: In Hollenstrain et al. (2013), they suggest that Eqn 49(a) in the
    // original Rubin and Attia paper should be used.

    // Method implemented below is translated from the SmoothMultiPhase fortran
    // subroutine of Rubin

    tmp<GeometricField<scalar, PatchField, Mesh>> tIbar
    (
        GeometricField<scalar, PatchField, Mesh>::New
        (
            "Ibar",
            mesh(),
            dimensionedScalar("zero", dimless, 1.0)
        )
    );

    GeometricField<scalar, PatchField, Mesh>& Ibar = tIbar.ref();

    // Calculate internal field
    forAll(Ibar, celli)
    {
        const scalar fac1 = 2.0/3.0*(devBEbar[celli] && devBEbar[celli]);
        if (mag(fac1) > small)
        {
            const scalar fac2 = 4.0*(1.0 - det(devBEbar[celli]))/pow(fac1, 1.5);

            if (fac2 >= 1.0)
            {
//                 scalar D = 1.0/sqr(fac1);
//                 Ibar[celli] =
//                     sqrt(fac1)
//                    *cosh
//                     (
//                         log((1.0 + sqrt(D))/sqrt(1.0 - D))/3.0
//                     );
                Ibar[celli] = sqrt(fac1)*cosh(acosh(fac2)/3.0);
            }
            else
            {
                Ibar[celli] = sqrt(fac1)*cos(acos(fac2)/3.0);
            }
        }
    }

    // Calculate boundary field
    typename GeometricField<scalar, PatchField, Mesh>::Boundary& bIbar =
        Ibar.boundaryFieldRef();
    forAll(Ibar.boundaryField(), patchi)
    {
        // Take reference to patch fields for efficiency
        scalarField& pIbar = bIbar[patchi];
        const symmTensorField& pdevBEbar = devBEbar.boundaryField()[patchi];

        forAll(pIbar, facei)
        {
            const scalar fac1 = 2.0/3.0*(pdevBEbar[facei] && pdevBEbar[facei]);
            if (mag(fac1) > small)
            {
                const scalar fac2 =
                    4.0*(1.0 - det(pdevBEbar[facei]))/pow(fac1, 1.5);

                if (fac2 >= 1.0)
                {
//                     scalar D = 1.0/sqr(fac1);
//                     Ibar[facei] =
//                         3.0*sqrt(fac1)
//                        *cosh
//                         (
//                             log((1.0 + sqrt(D))/sqrt(1.0 - D))/3.0
//                         );
                    Ibar[facei] = 3.0*sqrt(fac1)*cosh(acosh(fac2)/3.0);
                }
                else
                {
                    Ibar[facei] = sqrt(fac1)*cos(acos(fac2)/3.0);
                }
            }
        }
    }

    return tIbar;
}

// ************************************************************************* //
