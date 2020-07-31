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
{
    this->lookupAndInitialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFODETransportModel::~PDFODETransportModel()
{}


// * * * * * * * * * * * * Public Member functions * * * * * * * * * * * * * //

void Foam::PDFODETransportModel::clearODEFields()
{
    forAll(momentsOld_, stepi)
    {
        momentsOld_[stepi].clear();
        momentsOld_[stepi].resize(nMoments());
    }
    forAll(deltaMoments_, stepi)
    {
        deltaMoments_[stepi].clear();
        deltaMoments_[stepi].resize(nMoments());
    }
}

// ************************************************************************* //
