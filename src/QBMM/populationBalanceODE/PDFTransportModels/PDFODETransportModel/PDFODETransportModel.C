/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 Alberto Passalacqua
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

#include "PDFODETransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFODETransportModel::PDFODETransportModel
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    integrationSystem
    (
        IOobject::groupName("PDFODETransportModel", name),
        mesh
    ),
    name_(name),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFODETransportModel::~PDFODETransportModel()
{}


// * * * * * * * * * * * * Public Member functions * * * * * * * * * * * * * //

void Foam::PDFODETransportModel::initializeODEFields()
{
    momentsOld_ = PtrList<PtrList<volScalarField>>(nMoments());
    deltaMoments_ = PtrList<PtrList<volScalarField>>(nMoments());
    forAll(momentsOld_, mi)
    {
        momentsOld_.set(mi, new PtrList<volScalarField>());
        deltaMoments_.set(mi, new PtrList<volScalarField>());
    }
}


void Foam::PDFODETransportModel::clearODEFields()
{
    forAll(momentsOld_, mi)
    {
        this->clearOld(momentsOld_[mi]);
        this->clearDelta(deltaMoments_[mi]);
    }
}

// ************************************************************************* //
