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

#include "phaseFluxScheme.H"
#include "MUSCLReconstructionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
phaseFluxScheme::interpolateField
(
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const word& name
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;
    label nCmpts = pTraits<Type>::nComponents;

    bool isRho = f.member() == "rho";
    autoPtr<MUSCLReconstructionScheme<Type>> fLimiter
    (
        MUSCLReconstructionScheme<Type>::New(f, name)
    );

    tmp<fieldType> fOwnTmp(fLimiter->interpolateOwn());
    tmp<fieldType> fNeiTmp(fLimiter->interpolateNei());

    const fieldType& fOwn = fOwnTmp();
    const fieldType& fNei = fNeiTmp();

    tmp<fieldType> tmpf
    (
        new fieldType
        (
            IOobject
            (
                f.name() + "f",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensioned<Type>("0", f.dimensions(), pTraits<Type>::zero)
        )
    );
    fieldType& ff = tmpf.ref();

    forAll(fOwn, facei)
    {
        Type& fi = ff[facei];
        const Type& fiOwn = fOwn[facei];
        const Type& fiNei = fNei[facei];
        for (label i = 0; i < nCmpts; i++)
        {
            setComponent(fi, i) = interpolate
            (
                component(fiOwn, i),
                component(fiNei, i),
                isRho,
                facei
            );
        }
    }

    forAll(f.boundaryField(), patchi)
    {
        Field<Type>& pff = ff.boundaryFieldRef()[patchi];
        const Field<Type>& pfOwn = fOwn.boundaryField()[patchi];
        const Field<Type>& pfNei = fNei.boundaryField()[patchi];

        forAll(pff, facei)
        {
            Type& fi = pff[facei];
            const Type& fiOwn = pfOwn[facei];
            const Type& fiNei = pfNei[facei];
            for (label i = 0; i < nCmpts; i++)
            {
                setComponent(fi, i) =
                    interpolate
                    (
                        component(fiOwn, i),
                        component(fiNei, i),
                        isRho,
                        facei, patchi
                    );
            }
        }
    }
    return tmpf;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
