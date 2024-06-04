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
    Istream& is,
    const bool overwrite
)
:
    ReconstructionScheme<Type>(phi, is, overwrite),
    gradPhis_(0),
    bound_(true),
    extrapolate_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::LinearMUSCLReconstructionScheme<Type>::~LinearMUSCLReconstructionScheme()
{}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class Type>
void Foam::LinearMUSCLReconstructionScheme<Type>::constructGradPhis() const
{
    gradPhis_.setSize(pTraits<Type>::nComponents);
    tmp<fv::gradScheme<scalar>> lgradientScheme
    (
        fv::gradScheme<scalar>::New
        (
            this->mesh_,
            this->mesh_.schemes().grad
            (
                "limitedGrad(" + this->phi_.name() + ")"
            )
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


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::LinearMUSCLReconstructionScheme<Type>::interpolateOwn() const
{
    if (gradPhis_.size() != pTraits<Type>::nComponents)
    {
        constructGradPhis();
    }

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tphiOwn
    (
        GeometricField<Type, fvsPatchField, surfaceMesh>::New
        (
            this->ownName(),
            this->mesh_,
            dimensioned<Type>(this->phi_.dimensions(), Zero)
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& phiOwn = tphiOwn.ref();

    const labelList& owner = this->mesh_.owner();
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
        vector drOwn(fc[facei] - cc[own]);

        for (direction cmpti = 0; cmpti < pTraits<Type>::nComponents; cmpti++)
        {
            setComponent(phiOwn[facei], cmpti) =
                component(this->phi_[own], cmpti)
              + component(limOwn[facei], cmpti)
               *(drOwn & this->gradPhis_[cmpti][own]);
        }
    }
    if (bound_)
    {
        const labelList& neighbour = this->mesh_.neighbour();
        forAll(owner, facei)
        {
            label own = owner[facei];
            label nei = neighbour[facei];
            phiOwn[facei] =
                min
                (
                    phiOwn[facei],
                    max(this->phi_[own], (this->phi_[nei]))
                );
            phiOwn[facei] =
                max
                (
                    phiOwn[facei],
                    min(this->phi_[own], (this->phi_[nei]))
                );
        }
    }

    forAll(this->phi_.boundaryField(), patchi)
    {
        const fvPatch& patch = this->mesh_.boundary()[patchi];
        const fvPatchField<Type>& pphi = this->phi_.boundaryField()[patchi];
        Field<Type>& pphiOwn = phiOwn.boundaryFieldRef()[patchi];

        if (pphi.coupled())
        {
            Field<Type>& pphiOwn = phiOwn.boundaryFieldRef()[patchi];
            Field<Type> pphipOwn(pphi.patchInternalField());
            Field<Type> pphipNei(pphi.patchNeighbourField());

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
            if (bound_)
            {
                forAll(pphipOwn, facei)
                {
                    pphiOwn[facei] =
                        min
                        (
                            pphiOwn[facei],
                            max(pphipOwn[facei], pphipNei[facei])
                        );
                    pphiOwn[facei] =
                        max
                        (
                            pphiOwn[facei],
                            min(pphipOwn[facei], pphipNei[facei])
                        );
                }
            }
        }
        else if (!pphi.fixesValue() && extrapolate_)
        {
            Field<Type> pphiI(pphi.patchInternalField());
            vectorField pdelta(patch.fvPatch::delta());

            for
            (
                direction cmpti = 0;
                cmpti < pTraits<Type>::nComponents;
                cmpti++
            )
            {
                Field<vector> pgradPhi
                (
                    this->gradPhis_[cmpti].boundaryField()[patchi].patchInternalField()
                );

                forAll(pphi, facei)
                {
                    setComponent(pphiOwn[facei], cmpti) =
                        component(pphiI[facei], cmpti)
                      + (pdelta[facei] & pgradPhi[facei]);
                }
            }
        }
        else
        {
            pphiOwn = pphi;
        }
    }

    return tphiOwn;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::LinearMUSCLReconstructionScheme<Type>::interpolateNei() const
{
    if (gradPhis_.size() != pTraits<Type>::nComponents)
    {
        constructGradPhis();
    }

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tphiNei
    (
        GeometricField<Type, fvsPatchField, surfaceMesh>::New
        (
            this->neiName(),
            this->mesh_,
            dimensioned<Type>(this->phi_.dimensions(), Zero)
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& phiNei = tphiNei.ref();

    // const labelList& owner = this->mesh_.owner();
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
        label nei = neighbour[facei];
        vector drNei(fc[facei] - cc[nei]);

        for (direction cmpti = 0; cmpti < pTraits<Type>::nComponents; cmpti++)
        {
            setComponent(phiNei[facei], cmpti) =
                component(this->phi_[nei], cmpti)
              + component(limNei[facei], cmpti)
               *(drNei & this->gradPhis_[cmpti][nei]);
        }
    }
    if (bound_)
    {
        const labelList& owner = this->mesh_.owner();
        forAll(owner, facei)
        {
            label own = owner[facei];
            label nei = neighbour[facei];
            phiNei[facei] =
                min
                (
                    phiNei[facei],
                    max(this->phi_[own], (this->phi_[nei]))
                );
            phiNei[facei] =
                max
                (
                    phiNei[facei],
                    min(this->phi_[own], (this->phi_[nei]))
                );
        }
    }

    forAll(this->phi_.boundaryField(), patchi)
    {
        const fvPatch& patch = this->mesh_.boundary()[patchi];
        const fvPatchField<Type>& pphi = this->phi_.boundaryField()[patchi];
        Field<Type>& pphiNei = phiNei.boundaryFieldRef()[patchi];

        if (pphi.coupled())
        {
            Field<Type> pphipOwn(pphi.patchInternalField());
            Field<Type> pphipNei(pphi.patchNeighbourField());

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
            if (bound_)
            {
                forAll(pphipOwn, facei)
                {
                    pphiNei[facei] =
                        min
                        (
                            pphiNei[facei],
                            max(pphipOwn[facei], pphipNei[facei])
                        );
                    pphiNei[facei] =
                        max
                        (
                            pphiNei[facei],
                            min(pphipOwn[facei], pphipNei[facei])
                        );
                }
            }
        }
        else if (!pphi.fixesValue() && extrapolate_)
        {
            Field<Type> pphiI(pphi.patchInternalField());
            vectorField pdelta(patch.fvPatch::delta());

            for
            (
                direction cmpti = 0;
                cmpti < pTraits<Type>::nComponents;
                cmpti++
            )
            {
                Field<vector> pgradPhi
                (
                    this->gradPhis_[cmpti].boundaryField()[patchi].patchInternalField()
                );

                forAll(pphi, facei)
                {
                    setComponent(pphiNei[facei], cmpti) =
                        component(pphiI[facei], cmpti)
                      + (pdelta[facei] & pgradPhi[facei]);
                }
            }
        }
        else
        {
            pphiNei = pphi;
        }

    }

    return tphiNei;
}


// ************************************************************************* //
