/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2020
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "errorEstimator.H"
#include "coupledMaxErrorFvPatchScalarField.H"
#include "mappedWallFvPatch.H"
#include "mappedMovingWallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(errorEstimator, 0);
    defineRunTimeSelectionTable(errorEstimator, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::volScalarField& Foam::errorEstimator::lookupOrConstruct
(
    const fvMesh& mesh
) const
{
    if (!mesh.foundObject<volScalarField>("error"))
    {
        wordList boundaryTypes(mesh.boundaryMesh().size(), "zeroGradient");
        forAll(boundaryTypes, patchi)
        {
            if
            (
                isA<mappedWallFvPatch>(mesh.boundary()[patchi])
             || isA<mappedMovingWallFvPatch>(mesh.boundary()[patchi])
            )
            {
                boundaryTypes[patchi] = coupledMaxErrorFvPatchScalarField::typeName;
            }
            if (debug)
            {
                Pout<< "Patch:" << mesh.boundary()[patchi].patch().name() <<nl
                    << " cellType:" << boundaryTypes[patchi] << endl;
            }
        }

        volScalarField* fPtr
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("error", name_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                0.0,
                boundaryTypes
            )
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return mesh.lookupObjectRef<volScalarField>("error");
}

Foam::errorEstimator::errorEstimator
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
:
    mesh_(mesh),
    name_(name),
    error_(lookupOrConstruct(mesh)),
    lowerRefine_(0.0),
    lowerUnrefine_(0.0),
    upperRefine_(0.0),
    upperUnrefine_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimator::~errorEstimator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::errorEstimator::read(const dictionary& dict)
{
    lowerRefine_ = dict.lookup<scalar>("lowerRefineLevel");
    lowerUnrefine_ = dict.lookup<scalar>("unrefineLevel");
    upperRefine_ = dict.lookupOrDefault("upperRefineLevel", great);
    upperUnrefine_ = dict.lookupOrDefault("upperUnrefineLevel", great);
}


void Foam::errorEstimator::getFieldValue(const word& name, volScalarField& f) const
{

    this->getFieldValueType<scalar>(name, f);
    this->getFieldValueType<vector>(name, f);
    this->getFieldValueType<symmTensor>(name, f);
    this->getFieldValueType<sphericalTensor>(name, f);
    this->getFieldValueType<tensor>(name, f);
}


void Foam::errorEstimator::normalize(volScalarField& error)
{
    forAll(error, celli)
    {
        if
        (
            error[celli] < lowerUnrefine_
         || error[celli] > upperUnrefine_
        )
        {
            error[celli] = -1.0;
        }
        else if
        (
            error[celli] > lowerRefine_
         && error[celli] < upperRefine_
        )
        {
            error[celli] = 1.0;
        }
        else
        {
            error[celli] = 0.0;
        }
    }
}

// ************************************************************************* //
