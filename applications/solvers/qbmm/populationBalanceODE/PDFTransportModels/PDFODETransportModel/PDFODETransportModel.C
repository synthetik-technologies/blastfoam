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
    name_(name),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFODETransportModel::~PDFODETransportModel()
{}


// * * * * * * * * * * * * Public Member functions * * * * * * * * * * * * * //

void Foam::PDFODETransportModel::setODEFields
(
    const label nSteps,
    const boolList& storeFields,
    const boolList& storeDeltas
)
{
    oldIs_.resize(nSteps);
    deltaIs_.resize(nSteps);
    label fi = 0;
    for (label i = 0; i < nSteps; i++)
    {
        if (storeFields[i])
        {
            oldIs_[i] = fi;
            fi++;
        }
        else
        {
            oldIs_[i] = -1;
        }
    }
    nOld_ = fi;

    fi = 0;
    for (label i = 0; i < nSteps; i++)
    {
        if (storeDeltas[i])
        {
            deltaIs_[i] = fi;
            fi++;
        }
        else
        {
            deltaIs_[i] = -1;
        }
    }
    nDelta_ = fi;

    momentsOld_.resize(nOld_);
    deltaMoments_.resize(nDelta_);
    forAll(momentsOld_, stepi)
    {
        momentsOld_.set
        (
            stepi,
            new PtrList<volScalarField>(nMoments())
        );
    }
    forAll(deltaMoments_, stepi)
    {
        deltaMoments_.set
        (
            stepi,
            new PtrList<volScalarField>(nMoments())
        );
    }
}


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
