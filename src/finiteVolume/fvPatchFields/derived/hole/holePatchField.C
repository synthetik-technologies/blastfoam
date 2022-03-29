/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "holePatchField.H"
#include "noSlipFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(holePatchField, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::holePatchField::holePatchField
(
    const fvPatch& p
)
:
    patch_(p),
    LName_("undefined"),
    L_(p.size(), great),
    holeDensity_(p.size(), small),
    openFraction_(p.size(), 0),
    mask_(p.size(), 0.0)
{}


Foam::holePatchField::holePatchField
(
    const fvPatch& p,
    const dictionary& dict
)
:
    patch_(p),
    LName_(dict.lookup("LName")),
    L_("L", dict, p.size()),
    holeDensity_("holeDensity", dict, p.size()),
    openFraction_(p.size(), 0),
    mask_(p.size(), 0.0)
{}


Foam::holePatchField::holePatchField
(
    const holePatchField& hpf,
    const fvPatch& p,
    const fvPatchFieldMapper& mapper
)
:
    patch_(p),
    LName_(hpf.LName_),
    L_(mapper(hpf.L_)),
    holeDensity_(mapper(hpf.holeDensity_)),
    openFraction_(mapper(hpf.openFraction_)),
    mask_(mapper(hpf.mask_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::holePatchField::autoMap
(
    const fvPatchFieldMapper& m
)
{
     m(L_, L_);
     m(holeDensity_, holeDensity_);
     m(openFraction_, openFraction_);
     m(mask_, mask_);
}


void Foam::holePatchField::rmap
(
    const holePatchField& hpf,
    const labelList& addr
)
{
    L_.rmap(hpf.L_, addr);
    holeDensity_.rmap(hpf.holeDensity_, addr);
    openFraction_.rmap(hpf.openFraction_, addr);
    mask_.rmap(hpf.mask_, addr);
}


void Foam::holePatchField::updateCoeffs()
{
    const fvPatchField<scalar>& pl =
            patch_.lookupPatchField<volScalarField, scalar>(LName_);
    scalarField l
    (
        pl.patchInternalField()
//         pl.coupled()
//       ? max(pl.patchInternalField(), pl.patchInternalField())
//       : pl.patchInternalField()
    );
    mask_ = pos0(L_ - l);
    openFraction_ =
        min
        (
            Foam::constant::mathematical::pi*sqr(L_/2.0)*holeDensity_,
            1.0
        );
}


void Foam::holePatchField::write(Ostream& os) const
{
    writeEntry(os, "L", L_);
    writeEntry(os, "LName", LName_);
}


void Foam::holePatchField::removeBoundaryGradient
(
    const volVectorField& U,
    volVectorField& gradP
)
{
    volVectorField::Boundary& bgradP = gradP.boundaryFieldRef();
    forAll(U.boundaryField(), patchi)
    {
        if (isA<holePatchField>(U.boundaryField()[patchi]))
        {
            const polyPatch& patch = U.mesh().boundaryMesh()[patchi];
            const scalarField& mask =
                dynamicCast<const holePatchField>
                (
                    U.boundaryField()[patchi]
                ).mask();

            vectorField gradPI(bgradP[patchi].patchInternalField());
            const vectorField n(patch.faceNormals());

            const vectorField& pgradP(bgradP[patchi]);

            vectorField nGradP((pgradP & n)*n);
            vectorField tanGradp(pgradP - nGradP);
            bgradP[patchi] =
                (
                    mask*nGradP
                ) + tanGradp;
        }
    }
}

Foam::tmp<Foam::volVectorField>
Foam::holePatchField::removeBoundaryGradient
(
    const volVectorField& U,
    const tmp<volVectorField>& ctgradP
)
{
    // Transfer ownership
    tmp<volVectorField> tgradP = ctgradP;
    ctgradP.clear();

    volVectorField& gradP = tgradP.ref();
    volVectorField::Boundary& bgradP = gradP.boundaryFieldRef();
    forAll(U.boundaryField(), patchi)
    {
        if (isA<holePatchField>(U.boundaryField()[patchi]))
        {
            const polyPatch& patch = U.mesh().boundaryMesh()[patchi];
            const scalarField& mask =
                dynamicCast<const holePatchField>
                (
                    U.boundaryField()[patchi]
                ).mask();

            vectorField gradPI(bgradP[patchi].patchInternalField());
            const vectorField n(patch.faceNormals());

            const vectorField& pgradP(bgradP[patchi]);

            vectorField nGradP((pgradP & n)*n);
            vectorField tanGradp(pgradP - nGradP);
            bgradP[patchi] =
                (
                    mask*nGradP
                ) + tanGradp;
        }
    }
    return tgradP;
}


// ************************************************************************* //
