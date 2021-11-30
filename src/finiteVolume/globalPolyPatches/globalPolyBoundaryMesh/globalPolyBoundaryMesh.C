/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License

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

#include "globalPolyBoundaryMesh.H"
#include "coupledGlobalPolyPatch.H"
#include "pointMesh.H"
#include "IOdictionary.H"
#include "hashedWordList.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(globalPolyBoundaryMesh, 0);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::globalPolyBoundaryMesh::globalPolyBoundaryMesh
(
    const polyMesh& mesh
)
:
    GlobalPolyBoundaryMesh(mesh),
    interfacesDicts_()
{
    if (mesh.time().db().foundObject<IOdictionary>("regionProperties"))
    {
        interfacesDicts_ =
            HashTable<dictionary>
            (
                mesh.time().db().lookupObject<IOdictionary>
                (
                    "regionProperties"
                ).lookup("interfaces")
            );
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::globalPolyBoundaryMesh::~globalPolyBoundaryMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::globalPolyBoundaryMesh::movePoints()
{
    forAllIter
    (
        HashPtrTable<globalPolyPatch>,
        patches_,
        iter
    )
    {
        iter()->movePoints();
    }
    return true;
}


void Foam::globalPolyBoundaryMesh::updateMesh(const mapPolyMesh& mpm)
{
    forAllIter
    (
        HashPtrTable<globalPolyPatch>,
        patches_,
        iter
    )
    {
        iter()->updateMesh();
    }
}



void Foam::globalPolyBoundaryMesh::reorderPatches
(
    const labelUList& newToOld,
    const bool validBoundary
)
{}


void Foam::globalPolyBoundaryMesh::addPatch(const label patchi)
{}


void Foam::globalPolyBoundaryMesh::setDisplacementField
(
    const word& region,
    const word& name
)
{
    if (displacementFields_.found(region))
    {
        displacementFields_[region] = name;
    }
    else
    {
        displacementFields_.insert(region, name);
    }
}

// * * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * //

const Foam::globalPolyPatch&
Foam::globalPolyBoundaryMesh::operator[](const word& patchName) const
{
    return this->operator[](mesh().boundaryMesh()[patchName]);
}

const Foam::globalPolyPatch&
Foam::globalPolyBoundaryMesh::operator[](const polyPatch& pp) const
{
    if (!patches_.found(pp.name()))
    {
        dictionary dict;
        if (interfacesDicts_.found(mesh_.name()))
        {
            const dictionary& mDict = interfacesDicts_[mesh().name()];
            if (mDict.found(pp.name()))
            {
                dict = mDict.subDict(pp.name());
            }
        }
        if (displacementFields_.found(mesh_.name()))
        {
            if (!dict.found("displacementField"))
            {
                dict.add
                (
                    "displacementField",
                    displacementFields_[mesh_.name()]
                );
            }
        }
        patches_.insert
        (
            pp.name(),
            globalPolyPatch::New(dict, pp).ptr()
        );
    }

    return *patches_[pp.name()];
}


const Foam::globalPolyPatch&
Foam::globalPolyBoundaryMesh::operator[](const pointPatch& pp) const
{
    return this->operator[](mesh().boundaryMesh()[pp.name()]);
}


const Foam::coupledGlobalPolyPatch&
Foam::globalPolyBoundaryMesh::operator()(const word& patchName) const
{
    return this->operator()(mesh().boundaryMesh()[patchName]);
}

const Foam::coupledGlobalPolyPatch&
Foam::globalPolyBoundaryMesh::operator()(const polyPatch& pp) const
{
    if (patches_.found(pp.name()))
    {
        HashPtrTable<globalPolyPatch>::iterator iter = patches_.find(pp.name());
        if (!isA<coupledGlobalPolyPatch>(*iter()))
        {
            patches_.erase(iter);
        }
    }

    if (!patches_.found(pp.name()))
    {
        dictionary& dict =
            const_cast<dictionary&>
            (
                interfacesDicts_[mesh_.name()].subDict(pp.name())
            );
        if (displacementFields_.found(mesh_.name()))
        {
            if (!dict.found("displacementField"))
            {
                dict.add
                (
                    "displacementField",
                    displacementFields_[mesh_.name()]
                );
            }
        }
        patches_.insert
        (
            pp.name(),
            new coupledGlobalPolyPatch(dict, pp)
        );
    }


    return dynamicCast<const coupledGlobalPolyPatch>(*patches_[pp.name()]);
}


const Foam::coupledGlobalPolyPatch&
Foam::globalPolyBoundaryMesh::operator()(const pointPatch& pp) const
{
    return this->operator()(mesh().boundaryMesh()[pp.name()]);
}
// ************************************************************************* //
