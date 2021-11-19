/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "fluxSchemeBase.H"
#include "MUSCLReconstructionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fluxSchemeBase::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& f
) const
{
    return interpolateField<Type>(f, f.name());
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fluxSchemeBase::interpolate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fOwn,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fNei,
    const word& name
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;
    label nCmpts = pTraits<Type>::nComponents;

    tmp<fieldType> tmpf
    (
        fieldType::New
        (
            name + "f",
            mesh_,
            dimensioned<Type>("0", fOwn.dimensions(), pTraits<Type>::zero)
        )
    );
    fieldType& ff = tmpf.ref();
    const bool isDensity = (fOwn.dimensions() == dimDensity);

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
                isDensity,
                facei
            );
        }
    }

    forAll(fOwn.boundaryField(), patchi)
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
                        isDensity,
                        facei, patchi
                    );
            }
        }
    }
    return tmpf;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fluxSchemeBase::interpolateField
(
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const word& name
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    autoPtr<MUSCLReconstructionScheme<Type>> fLimiter
    (
        MUSCLReconstructionScheme<Type>::New(f, name)
    );

    tmp<fieldType> tfOwn(fLimiter->interpolateOwn());
    tmp<fieldType> tfNei(fLimiter->interpolateNei());

    return interpolate(tfOwn(), tfNei(), name);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
