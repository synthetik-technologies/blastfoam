/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "explicitNonLinearSolid.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * * //

void Foam::solidModels::explicitNonLinearSolid::updateWavespeeds()
{
    wavespeed_ =
        fvc::interpolate(sqrt(this->mechanical().elasticModulus()/this->rho()));
    sWavespeed_ =
        fvc::interpolate
        (
            sqrt(this->mechanical().shearModulus()/this->rho())
        );

    vector eigVal;
    tensor eigVec;

    surfaceTensorField Ff(this->Ff());
    surfaceTensorField Cf(Ff.T() & Ff);
    forAll(Cf, facei)
    {
        eigenStructure(Cf[facei], eigVal, eigVec);
        scalar s = sqrt(cmptMin(eigVal));
        this->wavespeed_[facei] /= s;
        this->sWavespeed_[facei] /= s;
    }

    const surfaceTensorField::Boundary& bCf(Cf.boundaryField());
    surfaceScalarField::Boundary& bwavespeed =
        this->wavespeed_.boundaryFieldRef();
    surfaceScalarField::Boundary& bsWavespeed =
        this->sWavespeed_.boundaryFieldRef();
    forAll(bCf, patchi)
    {
        forAll(bCf[patchi], facei)
        {
            eigenStructure(bCf[patchi][facei], eigVal, eigVec);
            scalar s = sqrt(cmptMin(eigVal));
            bwavespeed[patchi][facei] /= s;
            bsWavespeed[patchi][facei] /= s;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidModels::explicitNonLinearSolid::explicitNonLinearSolid
(
    const word& type,
    fvMesh& mesh,
    const nonLinearGeometry::nonLinearType nonLinear,
    const bool isSolid
)
:
    explicitSolid(type, mesh, nonLinear, isSolid)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
