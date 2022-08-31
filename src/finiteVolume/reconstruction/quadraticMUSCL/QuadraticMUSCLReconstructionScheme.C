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

#include "QuadraticMUSCLReconstructionScheme.H"
#include "gradScheme.H"


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::QuadraticMUSCLReconstructionScheme<Type>::QuadraticMUSCLReconstructionScheme
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    Istream& is
)
:
    ReconstructionScheme<Type>(phi, is),
    gradPhis_(pTraits<Type>::nComponents),
    hessPhis_(pTraits<Type>::nComponents)
{
    tmp<fv::gradScheme<scalar>> lgradientScheme
    (
        fv::gradScheme<scalar>::New
        (
            this->mesh_,
            this->mesh_.gradScheme("limitedGradMUSCL")
        )
    );
    tmp<fv::gradScheme<vector>> hgradientScheme
    (
        fv::gradScheme<vector>::New
        (
            this->mesh_,
            this->mesh_.gradScheme("limitedHessMUSCL")
        )
    );
    for (direction cmpti = 0; cmpti < pTraits<Type>::nComponents; cmpti++)
    {
        tmp<volScalarField> phiCmpt
        (
            volScalarField::New
            (
                this->phi_.name() + "_" + Foam::name(cmpti),
                this->phi_.component(cmpti)
            )
        );
        gradPhis_.set
        (
            cmpti,
            lgradientScheme().grad(phiCmpt)
        );
        hessPhis_.set
        (
            cmpti,
            hgradientScheme().grad(gradPhis_[cmpti])
        );
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::QuadraticMUSCLReconstructionScheme<Type>::~QuadraticMUSCLReconstructionScheme()
{}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::QuadraticMUSCLReconstructionScheme<Type>::interpolateOwn() const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> GeoField;
    tmp<GeoField> tphiOwn
    (
        GeoField::New
        (
            this->phi_.name() + "Own",
            this->mesh_,
            dimensioned<Type>(this->phi_.dimensions(), Zero)
        )
    );
    GeoField& phiOwn = tphiOwn.ref();

    const labelList& owner = this->mesh_.owner();
    const labelList& neighbour = this->mesh_.neighbour();
    const vectorField& cc = this->mesh_.C();
    const vectorField& fc = this->mesh_.Cf();

    tmp<GeoField> tlimOwn(this->calcLimiter(1.0));
    const GeoField& limOwn = tlimOwn();

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector drOwn(fc[facei] - cc[own]);

        for (direction cmpti = 0; cmpti < pTraits<Type>::nComponents; cmpti++)
        {
            setComponent(phiOwn[facei], cmpti) =
                component(this->phi_[own], cmpti)
              + component(limOwn[facei], cmpti)
               *(
                    (drOwn & this->gradPhis_[cmpti][own])
                  + ((drOwn & hessPhis_[cmpti][own]) & drOwn)*0.5
                );
        }

        // Hard limit to min/max of owner/neighbour values
        phiOwn[facei] =
            min
            (
                max(this->phi_[own], this->phi_[nei]),
                max
                (
                    phiOwn[facei],
                    min(this->phi_[own], this->phi_[nei])
                )
            );
    }

    typename GeoField::Boundary& bphiOwn = phiOwn.boundaryFieldRef();
    forAll(this->phi_.boundaryField(), patchi)
    {
        const fvPatch& patch = this->mesh_.boundary()[patchi];
        const fvPatchField<Type>& pphi = this->phi_.boundaryField()[patchi];
        if (patch.coupled())
        {
            Field<Type>& pphiOwn = bphiOwn[patchi];
            Field<Type> pphipOwn(pphi.patchInternalField());
            Field<Type> pphipNei(pphi.patchNeighbourField());

            const Field<Type>& plimOwn(limOwn.boundaryField()[patchi]);
            vectorField pdeltaOwn(patch.fvPatch::delta());

            for
            (
                direction cmpti = 0;
                cmpti < pTraits<Type>::nComponents;
                cmpti++
            )
            {
                Field<vector> pgradPhiOwn
                (
                    gradPhis_[cmpti].boundaryField()[patchi].patchInternalField()
                );
                Field<tensor> phessPhiOwn
                (
                    hessPhis_[cmpti].boundaryField()[patchi].patchInternalField()
                );

                forAll(pphipOwn, facei)
                {
                    setComponent(pphiOwn[facei], cmpti) =
                        component(pphipOwn[facei], cmpti)
                      + component(plimOwn[facei], cmpti)
                       *(
                            (pdeltaOwn[facei] & pgradPhiOwn[facei])
                          + (
                                (pdeltaOwn[facei] & phessPhiOwn[facei])
                              & pdeltaOwn[facei]
                            )*0.5
                        );
                }
            }

            // Hard limit to min/max of owner/neighbour values
            pphiOwn =
                min
                (
                    max(pphipOwn, pphipNei),
                    max
                    (
                        min(pphipOwn, pphipNei),
                        pphiOwn
                    )
                );
        }
        else
        {
            bphiOwn[patchi] = this->phi_.boundaryField()[patchi];
        }
    }

    return tphiOwn;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::QuadraticMUSCLReconstructionScheme<Type>::interpolateNei() const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> GeoField;
    tmp<GeoField> tphiNei
    (
        GeoField::New
        (
            this->phi_.name() + "Nei",
            this->mesh_,
            dimensioned<Type>(this->phi_.dimensions(), Zero)
        )
    );
    GeoField& phiNei = tphiNei.ref();

    const labelList& owner = this->mesh_.owner();
    const labelList& neighbour = this->mesh_.neighbour();
    const vectorField& cc = this->mesh_.C();
    const vectorField& fc = this->mesh_.Cf();

    tmp<GeoField> tlimNei
    (
        this->calcLimiter(-1.0)
    );
    const GeometricField<Type, fvsPatchField, surfaceMesh>& limNei = tlimNei();

    forAll(neighbour, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector drNei(fc[facei] - cc[nei]);
        for (direction cmpti = 0; cmpti < pTraits<Type>::nComponents; cmpti++)
        {
            setComponent(phiNei[facei], cmpti) =
                component(this->phi_[nei], cmpti)
              + component(limNei[facei], cmpti)
               *(
                    (drNei & gradPhis_[cmpti][nei])
                  + ((drNei & hessPhis_[cmpti][nei]) & drNei)*0.5
                );
        }

        // Hard limit to min/max of owner/neighbour values
        phiNei[facei] =
            min
            (
                max(this->phi_[own], this->phi_[nei]),
                max
                (
                    phiNei[facei],
                    min(this->phi_[own], this->phi_[nei])
                )
            );
    }

    typename GeoField::Boundary& bphiNei = phiNei.boundaryFieldRef();

    forAll(this->phi_.boundaryField(), patchi)
    {
        const fvPatch& patch = this->mesh_.boundary()[patchi];
        const fvPatchField<Type>& pphi = this->phi_.boundaryField()[patchi];
        if (patch.coupled())
        {
            Field<Type>& pphiNei = bphiNei[patchi];
            Field<Type> pphipOwn(pphi.patchInternalField());
            Field<Type> pphipNei(pphi.patchNeighbourField());

            const Field<Type>& plimNei(limNei.boundaryField()[patchi]);
            vectorField pdeltaNei(patch.fvPatch::delta() - patch.delta());

            for
            (
                direction cmpti = 0;
                cmpti < pTraits<Type>::nComponents;
                cmpti++
            )
            {
                Field<vector> pgradPhiNei
                (
                    gradPhis_[cmpti].boundaryField()[patchi].patchNeighbourField()
                );
                Field<tensor> phessPhiNei
                (
                    hessPhis_[cmpti].boundaryField()[patchi].patchNeighbourField()
                );

                forAll(pphipNei, facei)
                {
                    setComponent(pphiNei[facei], cmpti) =
                        component(pphipNei[facei], cmpti)
                      + component(plimNei[facei], cmpti)
                       *(
                            (pdeltaNei[facei] & pgradPhiNei[facei])
                          + (
                                (pdeltaNei[facei] & phessPhiNei[facei])
                              & pdeltaNei[facei]
                            )*0.5
                        );
                }
            }

            // Hard limit to min/max of owner/neighbour values
            pphiNei =
                min
                (
                    max(pphipOwn, pphipNei),
                    max
                    (
                        min(pphipOwn, pphipNei),
                        pphiNei
                    )
                );
        }
        else
        {
            bphiNei[patchi] = this->phi_.boundaryField()[patchi];
        }
    }

    return tphiNei;
}


// ************************************************************************* //
