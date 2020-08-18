/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
18-08-2020 Jeff Heylmun:    | Quadratic interpolation with MUSCL reconstruction
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

#include "quadraticMUSCLReconstructionScheme.H"
#include "gradScheme.H"


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::quadraticMUSCLReconstructionScheme<Type>::quadraticMUSCLReconstructionScheme
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    Istream& is
)
:
    MUSCLReconstructionScheme<Type>(phi, is),
    gradPhis_(pTraits<Type>::nComponents),
    hessPhis_(pTraits<Type>::nComponents)
{
    tmp<fv::gradScheme<scalar>> gradientScheme
    (
        fv::gradScheme<scalar>::New
        (
            this->mesh_,
            this->mesh_.gradScheme("grad(" + this->phi_.name() + ")")
        )
    );
    tmp<fv::gradScheme<scalar>> lgradientScheme
    (
        fv::gradScheme<scalar>::New
        (
            this->mesh_,
            this->mesh_.gradScheme("limitedGrad(" + this->phi_.name() + ")")
        )
    );
    tmp<fv::gradScheme<vector>> hgradientScheme
    (
        fv::gradScheme<vector>::New
        (
            this->mesh_,
            this->mesh_.gradScheme("limitedGrad(" + this->phi_.name() + ")")
        )
    );
    for (direction cmpti = 0; cmpti < pTraits<Type>::nComponents; cmpti++)
    {
        volScalarField phiCmpt
        (
            this->phi_.name() + "_" + Foam::name(cmpti),
            this->phi_.component(cmpti)
        );
        gradPhis_.set
        (
            cmpti,
            lgradientScheme().grad(phiCmpt)
        );
        hessPhis_.set
        (
            cmpti,
            hgradientScheme().grad(gradientScheme().grad(phiCmpt))
        );
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::quadraticMUSCLReconstructionScheme<Type>::~quadraticMUSCLReconstructionScheme()
{}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::quadraticMUSCLReconstructionScheme<Type>::interpolateOwn() const
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
    const labelList& neighbour = this->mesh_.neighbour();
    const vectorField& cc = this->mesh_.cellCentres();
    const vectorField& fc = this->mesh_.faceCentres();

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
               *(
                    (drOwn & this->gradPhis_[cmpti][own])
                  + ((drOwn & hessPhis_[cmpti][own]) & drOwn)
                );
        }

        // Hard limit to min/max of owner/neighbour values
        phiOwn[facei] = max(phiOwn[facei], minVal);
        phiOwn[facei] = min(phiOwn[facei], maxVal);
    }

    forAll(this->phi_.boundaryField(), patchi)
    {
        const fvPatch& patch = this->mesh_.boundary()[patchi];
        const fvPatchField<Type>& pphi = this->phi_.boundaryField()[patchi];
        if (patch.coupled())
        {
            Field<Type> pphiOwn(pphi.patchInternalField());
            Field<Type> pphiNei(pphi.patchNeighbourField());

            Field<Type> minVal(min(pphiOwn, pphiNei));
            Field<Type> maxVal(max(pphiOwn, pphiNei));

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
                    gradPhis_[cmpti].boundaryField()[patchi].patchInternalField()
                );
                Field<tensor> phessPhiOwn
                (
                    hessPhis_[cmpti].boundaryField()[patchi].patchInternalField()
                );

                forAll(pphiOwn, facei)
                {
                    setComponent(phiOwn.boundaryFieldRef()[patchi][facei], cmpti) =
                        component(pphiOwn[facei], cmpti)
                      + component(plimOwn[facei], cmpti)
                       *(
                            (pdeltaOwn[facei] & pgradPhiOwn[facei])
                          + (
                                (pdeltaOwn[facei] & phessPhiOwn[facei])
                              & pdeltaOwn[facei]
                            )
                        );
                }
            }

            // Hard limit to min/max of owner/neighbour values
            phiOwn.boundaryFieldRef()[patchi] =
                max(minVal, phiOwn.boundaryField()[patchi]);
            phiOwn.boundaryFieldRef()[patchi] =
                min(maxVal, phiOwn.boundaryField()[patchi]);
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
Foam::quadraticMUSCLReconstructionScheme<Type>::interpolateNei() const
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

    const labelList& owner = this->mesh_.owner();
    const labelList& neighbour = this->mesh_.neighbour();
    const vectorField& cc = this->mesh_.cellCentres();
    const vectorField& fc = this->mesh_.faceCentres();

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
               *(
                    (drNei & gradPhis_[cmpti][nei])
                  + ((drNei & hessPhis_[cmpti][nei]) & drNei)
                );
        }

        // Hard limit to min/max of owner/neighbour values
        phiNei[facei] = max(phiNei[facei], minVal);
        phiNei[facei] = min(phiNei[facei], maxVal);
    }

    forAll(this->phi_.boundaryField(), patchi)
    {
        const fvPatch& patch = this->mesh_.boundary()[patchi];
        const fvPatchField<Type>& pphi = this->phi_.boundaryField()[patchi];
        if (patch.coupled())
        {
            Field<Type> pphiOwn(pphi.patchInternalField());
            Field<Type> pphiNei(pphi.patchNeighbourField());

            Field<Type> minVal(min(pphiOwn, pphiNei));
            Field<Type> maxVal(max(pphiOwn, pphiNei));

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
                    gradPhis_[cmpti].boundaryField()[patchi].patchNeighbourField()
                );
                Field<tensor> phessPhiNei
                (
                    hessPhis_[cmpti].boundaryField()[patchi].patchNeighbourField()
                );

                forAll(pphiNei, facei)
                {
                    setComponent(phiNei.boundaryFieldRef()[patchi][facei], cmpti) =
                        component(pphiNei[facei], cmpti)
                      + component(plimNei[facei], cmpti)
                       *(
                            (pdeltaNei[facei] & pgradPhiNei[facei])
                          + (
                                (pdeltaNei[facei] & phessPhiNei[facei])
                              & pdeltaNei[facei]
                            )
                        );
                }
            }

            // Hard limit to min/max of owner/neighbour values
            phiNei.boundaryFieldRef()[patchi] =
                max(minVal, phiNei.boundaryField()[patchi]);
            phiNei.boundaryFieldRef()[patchi] =
                min(maxVal, phiNei.boundaryField()[patchi]);
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
