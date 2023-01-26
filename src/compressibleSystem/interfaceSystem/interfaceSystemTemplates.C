/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "interfaceSystem.H"

// * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::interfaceSystem::correctInterfaceField
(
    const word& schemeKey,
    const GeometricField<Type, fvPatchField, volMesh>& alpha,
    const GeometricField<Type, fvPatchField, volMesh>& psi,
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& nHatf,
    GeometricField<Type, fvsPatchField, surfaceMesh>& psif
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tpsiIntf
    (
        fvc::interpolate
        (
            psi,
            nHatf,
            reconstruction::scheme(schemeKey, psi.group(), psi.mesh())
        )
    );
    correctInterfaceField(alpha, tpsiIntf(), psif);
    return tpsiIntf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::interfaceSystem::correctInterfaceField
(
    const word& schemeKey,
    const GeometricField<Type, fvPatchField, volMesh>& alpha,
    const GeometricField<Type, fvPatchField, volMesh>& psi,
    const tmp<GeometricField<scalar, fvsPatchField, surfaceMesh>>& nHatf,
    GeometricField<Type, fvsPatchField, surfaceMesh>& psif
)
{
    return correctInterfaceField(schemeKey, alpha, psi, nHatf(), psif);
}


template<class Type>
void Foam::interfaceSystem::correctInterfaceField
(
    const GeometricField<Type, fvPatchField, volMesh>& alpha,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& psiIntf,
    GeometricField<Type, fvsPatchField, surfaceMesh>& psif
)
{
    const labelList& owner = alpha.mesh().owner();
    const labelList& neighbour = alpha.mesh().neighbour();
    forAll(psif, facei)
    {
        // Only apply to faces where
        const label own = owner[facei];
        const label nei = neighbour[facei];
        if (0.99 > alpha[own] || 0.99 > alpha[nei])
        {
            psif[facei] = psiIntf[facei];
        }
    }

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary& bpsif =
        psif.boundaryFieldRef();
    forAll(bpsif, patchi)
    {
        fvsPatchField<Type>& ppsif = bpsif[patchi];
        const fvsPatchField<Type>& ppsiIntf = psiIntf.boundaryField()[patchi];
        const fvPatchField<scalar>& palpha = alpha.boundaryField()[patchi];
        if (ppsif.coupled())
        {
            Field<Type> alphaOwn(palpha.patchInternalField());
            Field<Type> alphaNei(palpha.patchNeighbourField());
            forAll(ppsif, facei)
            {
                // Only apply to faces where
                if (0.99 > alphaOwn[facei] || 0.99 > alphaNei[facei])
                {
                    ppsif[facei] = ppsiIntf[facei];
                }
            }
        }
    }
}


template<class Type>
void Foam::interfaceSystem::correctInterfaceField
(
    const GeometricField<Type, fvPatchField, volMesh>& alpha,
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& psiIntf,
    GeometricField<Type, fvsPatchField, surfaceMesh>& psif
)
{
    correctInterfaceField(alpha, psiIntf(), psif);
    psiIntf.clear();
}



template<class Type>
void Foam::interfaceSystem::correctInterfaceField
(
    const word& schemeKey,
    const GeometricField<Type, fvPatchField, volMesh>& alpha,
    const GeometricField<Type, fvPatchField, volMesh>& psi,
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& nHatf,
    GeometricField<Type, fvsPatchField, surfaceMesh>& psiOwn,
    GeometricField<Type, fvsPatchField, surfaceMesh>& psiNei
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tpsiIntf
    (
        correctInterfaceField
        (
            schemeKey,
            alpha,
            psi,
            nHatf,
            psiOwn
        )
    );
    correctInterfaceField(alpha, tpsiIntf, psiNei);
}
// ************************************************************************* //
