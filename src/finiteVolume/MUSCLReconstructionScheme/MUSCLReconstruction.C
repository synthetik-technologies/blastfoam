/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
18-08-2020 Jeff Heylmun:    | MUSCL reconstruction
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

#include "MUSCLReconstruction.H"
#include "gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template
<
    class Type,
    class MUSCLType,
    class Limiter,
    template<class> class LimitFunc
>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::MUSCLReconstruction<Type, MUSCLType, Limiter, LimitFunc>::calcLimiter
(
    const scalar& dir
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;
    const word limiterFieldName(type() + "Limiter(" + this->phi_.name() + ')');
    tmp<fieldType> tlimiterField
    (
        new fieldType
        (
            IOobject
            (
                limiterFieldName,
                this->mesh_.time().timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimensioned<Type>(dimless, Zero)
        )
    );
    fieldType& limiterField = tlimiterField.ref();

    const surfaceScalarField& CDweights =
        this->mesh_.surfaceInterpolation::weights();

    const labelUList& owner = this->mesh_.owner();
    const labelUList& neighbour = this->mesh_.neighbour();

    const vectorField& C = this->mesh_.C();

    tmp<fv::gradScheme<scalar>> gradientScheme
    (
        fv::gradScheme<scalar>::New
        (
            this->mesh_,
            this->mesh_.gradScheme(word("grad(" + this->phi_.name() + ")"))
        )
    );

    for (direction cmpti = 0; cmpti < pTraits<Type>::nComponents; cmpti++)
    {
        volScalarField phiCmpt(this->phi_.component(cmpti));
        tmp<GeometricField<typename Limiter::phiType, fvPatchField, volMesh>>
            tlPhi = LimitFunc<scalar>()(phiCmpt);

        const GeometricField<typename Limiter::phiType, fvPatchField, volMesh>&
            lPhi = tlPhi();

        tmp<GeometricField<typename Limiter::gradPhiType, fvPatchField, volMesh>>
            tgradc(gradientScheme().grad(lPhi));
        const GeometricField<typename Limiter::gradPhiType, fvPatchField, volMesh>&
            gradc = tgradc();

        forAll(owner, face)
        {
            label own = owner[face];
            label nei = neighbour[face];


            setComponent(limiterField[face], cmpti) =
                Limiter::limiter
                (
                    CDweights[face],
                    dir,
                    lPhi[own],
                    lPhi[nei],
                    gradc[own],
                    gradc[nei],
                    C[nei] - C[own]
                );
        }

        typename fieldType::Boundary& bLim = limiterField.boundaryFieldRef();

        forAll(bLim, patchi)
        {
            Field<Type>& pLim = bLim[patchi];

            if (bLim[patchi].coupled())
            {
                const scalarField& pCDweights = CDweights.boundaryField()[patchi];

                const Field<typename Limiter::phiType> plPhiP
                (
                    lPhi.boundaryField()[patchi].patchInternalField()
                );
                const Field<typename Limiter::phiType> plPhiN
                (
                    lPhi.boundaryField()[patchi].patchNeighbourField()
                );
                const Field<typename Limiter::gradPhiType> pGradcP
                (
                    gradc.boundaryField()[patchi].patchInternalField()
                );
                const Field<typename Limiter::gradPhiType> pGradcN
                (
                    gradc.boundaryField()[patchi].patchNeighbourField()
                );

                // Build the d-vectors
                vectorField pd
                (
                    CDweights.boundaryField()[patchi].patch().delta()
                );

                forAll(pLim, facei)
                {
                    setComponent(pLim[facei], cmpti) =
                        Limiter::limiter
                        (
                            pCDweights[facei],
                            dir,
                            plPhiP[facei],
                            plPhiN[facei],
                            pGradcP[facei],
                            pGradcN[facei],
                            pd[facei]
                        );
                }
            }
            else
            {
                forAll(pLim, facei)
                {
                    setComponent(pLim[facei], cmpti) = 1.0;
                }
            }
        }
    }

    return tlimiterField;
}



// ************************************************************************* //
