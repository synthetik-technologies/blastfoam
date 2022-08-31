/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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

#include "speciesMixtureField.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::speciesMixtureField<ThermoType>::speciesMixtureField
(
    const fvMesh& mesh,
    const PtrList<volScalarField>& Ys,
    const PtrList<ThermoType>& speciesData,
    const word& phaseName
)
:
    PtrList<ThermoType>(mesh.nCells()),
    MeshObject
    <
        fvMesh,
        DistributeableMeshObject,
        speciesMixtureField<ThermoType>
    >
    (
        mesh,
        IOobject
        (
            IOobject::groupName("speciesThermo", phaseName),
            mesh.time().timeName(),
            mesh
        )
    ),
    mesh_(mesh),
    Ys_(Ys),
    speciesData_(speciesData),
    mixture_(speciesData[0]),
    faceMixtures_(nBoundaryFaces()),
    start_(Ys_[0].boundaryField().size())
{
    label i = 0;
    forAll(start_, patchi)
    {
        start_[patchi] = i;
        i += Ys_[0].boundaryField()[patchi].size();
    }

    //- Allocate cell mixtures
    forAll(*this, celli)
    {
        this->set
        (
            celli,
            new ThermoType(cellMixture(celli))
        );
    }

    //- Allocate boundary mixtures
    forAll(Ys_[0].boundaryField(), patchi)
    {
        forAll(Ys_[0].boundaryField()[patchi], facei)
        {
            faceMixtures_.set
            (
                patchFaceIndex(patchi, facei),
                new ThermoType(patchFaceMixture(patchi, facei))
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::speciesMixtureField<ThermoType>::~speciesMixtureField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::speciesMixtureField<ThermoType>::updateCells()
{
    //- Resize and update cell mixtures
    const label nCellsOld = this->size();
    const label nCells = mesh_.nCells();
    const label nC(min(nCells, nCellsOld));
    this->resize(nCells);

    for (label celli = 0; celli < nC; celli++)
    {
        this->operator[](celli) = cellMixture(celli);
    }
    for (label celli = nCellsOld; celli < nCells; celli++)
    {
        this->set(celli, new ThermoType(cellMixture(celli)));
    }
}


template<class ThermoType>
void Foam::speciesMixtureField<ThermoType>::updatePatches
(
    const bool validBoundary
)
{
    label i = 0;
    const fvBoundaryMesh& bm = mesh_.boundary();
    start_.resize(bm.size());
    forAll(bm, patchi)
    {
        start_[patchi] = i;
        i += bm[patchi].size();
    }

    //- Resize boundary
    faceMixtures_.resize(nBoundaryFaces());

    //- Update boundaries that were previously allocated
    forAll(bm, patchi)
    {
        forAll(bm[patchi], facei)
        {
            label pfi = patchFaceIndex(patchi, facei);
            if (faceMixtures_.set(pfi))
            {
                faceMixtures_[pfi] =
                    validBoundary
                  ? patchFaceMixture(patchi, facei)
                  : mixture_;
            }
            else
            {
                faceMixtures_.set
                (
                    pfi,
                    validBoundary
                  ? new ThermoType(patchFaceMixture(patchi, facei))
                  : new ThermoType(mixture_)
                );
            }
        }
    }
}

template<class ThermoType>
void Foam::speciesMixtureField<ThermoType>::updateMixture()
{
    PtrList<ThermoType>& cells(*this);
    forAll(cells, celli)
    {
        cells[celli] = cellMixture(celli);
    }

    forAll(Ys_[0].boundaryField(), patchi)
    {
        forAll(Ys_[0].boundaryField()[patchi], facei)
        {
            faceMixtures_[patchFaceIndex(patchi, facei)] =
                patchFaceMixture(patchi, facei);
        }
    }
}


template<class ThermoType>
const ThermoType& Foam::speciesMixtureField<ThermoType>::cellMixture
(
    const label celli
) const
{
    mixture_ = Ys_[0][celli]*speciesData_[0];

    for (label n = 1; n < Ys_.size(); n++)
    {
        mixture_ += Ys_[n][celli]*speciesData_[n];
    }

    return mixture_;
}


template<class ThermoType>
const ThermoType&
Foam::speciesMixtureField<ThermoType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    mixture_ = Ys_[0].boundaryField()[patchi][facei]*speciesData_[0];

    for (label n = 1; n < Ys_.size(); n++)
    {
        mixture_ += Ys_[n].boundaryField()[patchi][facei]*speciesData_[n];
    }

    return mixture_;
}


// ************************************************************************* //
