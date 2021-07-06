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
#include "timeControlFunctionObject.H"
#include "probes.H"
#include "blastProbes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(errorEstimator, 0);
    defineRunTimeSelectionTable(errorEstimator, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::errorEstimator::constructError
(
    const fvMesh& mesh
) const
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

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("error", name_),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                debug ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
            ),
            mesh,
            0.0,
            boundaryTypes
        )
    );
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
    error_(constructError(mesh)),
    lowerRefine_(0.0),
    lowerUnrefine_(0.0),
    upperRefine_(0.0),
    upperUnrefine_(0.0),
    maxLevel_(10),
    refineProbes_(dict.lookupOrDefault("refineProbes", true))
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
    maxLevel_ = dict.lookup<label>("maxRefinement");
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

    if (!refineProbes_)
    {
        return;
    }

    const functionObjectList& funcs(mesh_.time().functionObjects());
    labelList map;
    forAll(funcs, i)
    {
        if (isA<functionObjects::timeControl>(funcs[i]))
        {
            const functionObjects::timeControl& tc
            (
                dynamicCast<const functionObjects::timeControl>(funcs[i])
            );
            bool sameRegion = true;
            const dictionary& dict = tc.dict();
            if (dict.found("region"))
            {
                word regionName = dict.lookup<word>("region");
                sameRegion = regionName == error_.mesh().name();
            }

            if (isA<probes>(tc.filter()) && sameRegion)
            {
                const probes& p(refCast<const probes>(tc.filter()));
                const labelList& cells(p.elements());
                forAll(cells, j)
                {
                    label celli = cells[j];
                    if (celli >= 0)
                    {
                        error[celli] = 1.0;
                    }
                }
            }
            if (isA<blastProbes>(tc.filter()) && sameRegion)
            {
                const blastProbes& p(refCast<const blastProbes>(tc.filter()));
                const labelList& cells(p.elements());
                forAll(cells, j)
                {
                    label celli = cells[j];
                    if (celli >= 0)
                    {
                        error[celli] = 1.0;
                    }
                }
            }
        }
    }
}

Foam::labelList Foam::errorEstimator::maxRefinement() const
{
    return labelList(mesh_.nCells(), maxLevel_);
}


// ************************************************************************* //
