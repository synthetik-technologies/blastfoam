/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
18-08-2020 Jeff Heylmun:    | Upwind MUSCL reconstruction
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "upwindMUSCLReconstructionScheme.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::upwindMUSCLReconstructionScheme<Type>::upwindMUSCLReconstructionScheme
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    Istream& is
)
:
    MUSCLReconstructionScheme<Type>(phi, is)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::upwindMUSCLReconstructionScheme<Type>::~upwindMUSCLReconstructionScheme()
{}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::upwindMUSCLReconstructionScheme<Type>::interpolateOwn() const
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tphiOwn
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                this->phi_.name() + "Own",
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensioned<Type>(this->phi_.dimensions(), Zero)
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& phiOwn = tphiOwn.ref();

    const labelList& owner = this->mesh_.owner();
    forAll(owner, facei)
    {
        phiOwn[facei] = this->phi_[owner[facei]];
    }

    forAll(this->phi_.boundaryField(), patchi)
    {
        const fvPatch& patch = this->mesh_.boundary()[patchi];
        const fvPatchField<Type>& pphi = this->phi_.boundaryField()[patchi];
        if (patch.coupled())
        {
            phiOwn.boundaryFieldRef()[patchi] = pphi.patchInternalField();
        }
        else
        {
            phiOwn.boundaryFieldRef()[patchi] =
                this->phi_.boundaryField()[patchi];
        }
    }

    return tphiOwn;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::upwindMUSCLReconstructionScheme<Type>::interpolateNei() const
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tphiNei
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                this->phi_.name() + "Nei",
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensioned<Type>(this->phi_.dimensions(), Zero)
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& phiNei = tphiNei.ref();

    const labelList& nei = this->mesh_.neighbour();
    forAll(nei, facei)
    {
        phiNei[facei] = this->phi_[nei[facei]];
    }

    forAll(this->phi_.boundaryField(), patchi)
    {
        const fvPatch& patch = this->mesh_.boundary()[patchi];
        const fvPatchField<Type>& pphi = this->phi_.boundaryField()[patchi];
        if (patch.coupled())
        {
            phiNei.boundaryFieldRef()[patchi] = pphi.patchNeighbourField();
        }
        else
        {
            phiNei.boundaryFieldRef()[patchi] =
                this->phi_.boundaryField()[patchi];
        }
    }

    return tphiNei;
}

// ************************************************************************* //
