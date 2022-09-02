/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020
     \\/     M anipulation  | Synthetik Applied Technology
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

#include "LinearMUSCLReconstructionScheme.H"
#include "gradScheme.H"


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::LinearMUSCLReconstructionScheme<Type>::LinearMUSCLReconstructionScheme
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    Istream& is
)
:
    ReconstructionScheme<Type>(phi, is),
    gradPhis_(pTraits<Type>::nComponents)
{
    tmp<fv::gradScheme<scalar>> lgradientScheme
    (
        fv::gradScheme<scalar>::New
        (
            this->mesh_,
            this->mesh_.gradScheme("limitedGrad(" + this->phi_.name() + ")")
        )
    );
    for (direction cmpti = 0; cmpti < pTraits<Type>::nComponents; cmpti++)
    {
        gradPhis_.set
        (
            cmpti,
            lgradientScheme().grad(this->phi_.component(cmpti))
        );
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::LinearMUSCLReconstructionScheme<Type>::~LinearMUSCLReconstructionScheme()
{}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::LinearMUSCLReconstructionScheme<Type>::interpolateOwn() const
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tphiOwn
    (
        GeometricField<Type, fvsPatchField, surfaceMesh>::New
        (
            this->phi_.name() + "Own",
            this->mesh_,
            dimensioned<Type>(this->phi_.dimensions(), Zero)
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& phiOwn = tphiOwn.ref();

    const labelList& owner = this->mesh_.owner();
    const labelList& neighbour = this->mesh_.neighbour();
    const vectorField& cc = this->mesh_.C();
    const vectorField& fc = this->mesh_.Cf();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tlimOwn
    (
        this->calcLimiter(1.0)
    );
    const GeometricField<Type, fvsPatchField, surfaceMesh>& limOwn = tlimOwn();

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        Type minVal(min(this->phi_[own], this->phi_[nei]));
        Type maxVal(max(this->phi_[own], this->phi_[nei]));

        vector drOwn(fc[facei] - cc[own]);

        for (direction cmpti = 0; cmpti < pTraits<Type>::nComponents; cmpti++)
        {
            setComponent(phiOwn[facei], cmpti) =
                component(this->phi_[own], cmpti)
              + component(limOwn[facei], cmpti)
               *(drOwn & this->gradPhis_[cmpti][own]);
        }

        // Hard limit to min/max of owner/neighbour values
        phiOwn[facei] = max(phiOwn[facei], minVal);
        phiOwn[facei] = min(phiOwn[facei], maxVal);
    }

    forAll(this->phi_.boundaryField(), patchi)
    {
        const fvPatch& patch = this->mesh_.boundary()[patchi];
        const fvPatchField<Type>& pphi = this->phi_.boundaryField()[patchi];
        if (pphi.coupled())
        {
            Field<Type>& pphiOwn = phiOwn.boundaryFieldRef()[patchi];
            Field<Type> pphipOwn(pphi.patchInternalField());
            Field<Type> pphipNei(pphi.patchNeighbourField());

            Field<Type> minVal(min(pphipOwn, pphipNei));
            Field<Type> maxVal(max(pphipOwn, pphipNei));

            const Field<Type>& plimOwn
            (
                limOwn.boundaryField()[patchi]
            );

            vectorField pdeltaOwn
            (
                patch.fvPatch::delta()
            );

            for
            (
                direction cmpti = 0;
                cmpti < pTraits<Type>::nComponents;
                cmpti++
            )
            {
                Field<vector> pgradPhiOwn
                (
                    this->gradPhis_[cmpti].boundaryField()[patchi].patchInternalField()
                );

                forAll(pphipOwn, facei)
                {
                    setComponent(pphiOwn[facei], cmpti) =
                        component(pphipOwn[facei], cmpti)
                      + component(plimOwn[facei], cmpti)
                       *(pdeltaOwn[facei] & pgradPhiOwn[facei]);
                }
            }

            // Hard limit to min/max of owner/neighbour values
            pphiOwn = max(minVal, phiOwn.boundaryField()[patchi]);
            pphiOwn = min(maxVal, phiOwn.boundaryField()[patchi]);
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
Foam::LinearMUSCLReconstructionScheme<Type>::interpolateNei() const
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tphiNei
    (
        GeometricField<Type, fvsPatchField, surfaceMesh>::New
        (
            this->phi_.name() + "Nei",
            this->mesh_,
            dimensioned<Type>(this->phi_.dimensions(), Zero)
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& phiNei = tphiNei.ref();

    const labelList& owner = this->mesh_.owner();
    const labelList& neighbour = this->mesh_.neighbour();
    const vectorField& cc = this->mesh_.C();
    const vectorField& fc = this->mesh_.Cf();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tlimNei
    (
        this->calcLimiter(-1.0)
    );
    const GeometricField<Type, fvsPatchField, surfaceMesh>& limNei = tlimNei();

    forAll(neighbour, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        Type minVal(min(this->phi_[own], this->phi_[nei]));
        Type maxVal(max(this->phi_[own], this->phi_[nei]));

        vector drNei(fc[facei] - cc[nei]);
        for (direction cmpti = 0; cmpti < pTraits<Type>::nComponents; cmpti++)
        {
            setComponent(phiNei[facei], cmpti) =
                component(this->phi_[nei], cmpti)
              + component(limNei[facei], cmpti)
               *(drNei & this->gradPhis_[cmpti][nei]);
        }

        // Hard limit to min/max of owner/neighbour values
        phiNei[facei] = max(phiNei[facei], minVal);
        phiNei[facei] = min(phiNei[facei], maxVal);
    }

    forAll(this->phi_.boundaryField(), patchi)
    {
        const fvPatch& patch = this->mesh_.boundary()[patchi];
        const fvPatchField<Type>& pphi = this->phi_.boundaryField()[patchi];
        if (pphi.coupled())
        {
            Field<Type>& pphiNei = phiNei.boundaryFieldRef()[patchi];
            Field<Type> pphipOwn(pphi.patchInternalField());
            Field<Type> pphipNei(pphi.patchNeighbourField());

            Field<Type> minVal(min(pphipOwn, pphipNei));
            Field<Type> maxVal(max(pphipOwn, pphipNei));

            const Field<Type>& plimNei
            (
                limNei.boundaryField()[patchi]
            );
            vectorField pdeltaNei
            (
                patch.fvPatch::delta() - patch.delta()
            );

            for
            (
                direction cmpti = 0;
                cmpti < pTraits<Type>::nComponents;
                cmpti++
            )
            {
                Field<vector> pgradPhiNei
                (
                    this->gradPhis_[cmpti].boundaryField()[patchi].patchNeighbourField()
                );

                forAll(pphipNei, facei)
                {
                    setComponent(pphiNei[facei], cmpti) =
                        component(pphipNei[facei], cmpti)
                      + component(plimNei[facei], cmpti)
                       *(pdeltaNei[facei] & pgradPhiNei[facei]);
                }
            }

            // Hard limit to min/max of owner/neighbour values
            pphiNei = max(minVal, phiNei.boundaryField()[patchi]);
            pphiNei = min(maxVal, phiNei.boundaryField()[patchi]);
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
