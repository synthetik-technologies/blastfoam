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
    UpdateableMeshObject<polyMesh>
    (
        IOobject::groupName("speciesThermo", phaseName),
        mesh
    ),
    mesh_(mesh),
    Ys_(Ys),
    speciesData_(speciesData),
    mixture_(speciesData[0]),
    faceMixtures_(mesh.boundary().size())
{
    forAll(*this, celli)
    {
        this->set
        (
            celli,
            new ThermoType(cellMixture(celli))
        );
    }
    forAll(faceMixtures_, patchi)
    {
        faceMixtures_.set
        (
            patchi,
            new PtrList<ThermoType>(mesh.boundary()[patchi].size())
        );
        forAll(faceMixtures_[patchi], facei)
        {
            faceMixtures_[patchi].set
            (
                facei,
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
void Foam::speciesMixtureField<ThermoType>::updateMesh(const mapPolyMesh& mpm)
{
    //- Only allocate list, updating is done later
    label nOld = this->size();
    label n = mesh_.nCells();
    this->resize(n);

    for (label celli = nOld; celli < n; celli++)
    {
        this->set
        (
            celli,
            new ThermoType(mixture_)
        );
    }

    nOld = faceMixtures_.size();
    n = mesh_.boundary().size();
    faceMixtures_.resize(n);
    for (label patchi = 0; patchi < min(nOld, n); patchi++)
    {
        label nFaces = mesh_.boundary()[patchi].size();
        label nFacesOld = faceMixtures_[patchi].size();
        faceMixtures_[patchi].resize(nFaces);
        for (label facei = nFacesOld; facei < nFaces; facei++)
        {
            faceMixtures_[patchi].set
            (
                facei,
                new ThermoType(mixture_)
            );
        }
    }
    for (label patchi = nOld; patchi < n; patchi++)
    {
        label nFaces = mesh_.boundary()[patchi].size();
        faceMixtures_.set
        (
            patchi,
            new PtrList<ThermoType>(nFaces)
        );
        forAll(faceMixtures_[patchi], facei)
        {
            faceMixtures_[patchi].set
            (
                facei,
                new ThermoType(mixture_)
            );
        }
    }
}

template<class ThermoType>
void Foam::speciesMixtureField<ThermoType>::updateCellMixtures()
{
    PtrList<ThermoType>& cells(*this);
    forAll(cells, celli)
    {
        cells[celli] = cellMixture(celli);
    }
}


template<class ThermoType>
void Foam::speciesMixtureField<ThermoType>::updateFaceMixtures
(
    const label patchi
)
{
    forAll(faceMixtures_[patchi], facei)
    {
        faceMixtures_[patchi][facei] = patchFaceMixture(patchi, facei);
    }
}



template<class ThermoType>
void Foam::speciesMixtureField<ThermoType>::update()
{
    updateCellMixtures();

    forAll(faceMixtures_, patchi)
    {
        updateFaceMixtures(patchi);
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
