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
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(globalPolyBoundaryMesh, 0);
}

bool Foam::globalPolyBoundaryMesh::clearOnMovement = true;

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::globalPolyBoundaryMesh::globalPolyBoundaryMesh
(
    const polyMesh& mesh
)
:
    GlobalPolyBoundaryMesh(mesh),
    interfaceDicts_(),
    readFromRP_(false)
{
    if (mesh.time().db().foundObject<IOdictionary>("regionProperties"))
    {
        const IOdictionary& regionProperties =
            mesh.time().db().lookupObject<IOdictionary>("regionProperties");
        if (regionProperties.found("interfaces"))
        {
            interfaceDicts_ =
                HashTable<dictionary>(regionProperties.lookup("interfaces"));
            readFromRP_ = true;
        }
    }
}


Foam::globalPolyBoundaryMesh::globalPolyBoundaryMesh
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    GlobalPolyBoundaryMesh(mesh),
    interfaceDicts_(dict.lookupOrDefault("interfaces", HashTable<dictionary>())),
    readFromRP_(false)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::globalPolyBoundaryMesh::~globalPolyBoundaryMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::globalPolyBoundaryMesh::isGlobal(const polyPatch& pp) const
{
    return patches_.found(pp.name());
}


bool Foam::globalPolyBoundaryMesh::isCoupled(const polyPatch& pp) const
{
    const polyMesh& mesh = pp.boundaryMesh().mesh();
    if (!interfaceDicts_.found(mesh.name()))
    {
        return false;
    }
    return interfaceDicts_[mesh.name()].found(pp.name());
}


Foam::labelList Foam::globalPolyBoundaryMesh::coupledPatches() const
{
    if (!interfaceDicts_.found(this->mesh().name()))
    {
        return labelList();
    }
    const dictionary& dict = interfaceDicts_[this->mesh().name()];
    DynamicList<label> patches(this->mesh().boundaryMesh().size());
    forAll(this->mesh().boundaryMesh(), patchi)
    {
        if (dict.isDict(this->mesh().boundaryMesh()[patchi].name()))
        {
            patches.append(this->mesh().boundaryMesh()[patchi].index());
        }
    }
    return patches;
}


void Foam::globalPolyBoundaryMesh::update()
{
    forAllIter
    (
        HashPtrTable<globalPolyPatch>,
        patches_,
        iter
    )
    {
        iter()->update();
    }
}


bool Foam::globalPolyBoundaryMesh::movePoints()
{
    forAllIter
    (
        HashPtrTable<globalPolyPatch>,
        patches_,
        iter
    )
    {
        iter()->movePoints(clearOnMovement);
    }
    return true;
}


void Foam::globalPolyBoundaryMesh::distribute(const polyDistributionMap& map)
{
    forAllIter
    (
        HashPtrTable<globalPolyPatch>,
        patches_,
        iter
    )
    {
        iter()->distribute(map);
    }
}


void Foam::globalPolyBoundaryMesh::topoChange(const polyTopoChangeMap& map)
{
    forAllIter
    (
        HashPtrTable<globalPolyPatch>,
        patches_,
        iter
    )
    {
        iter()->topoChange(map);
    }
}


void Foam::globalPolyBoundaryMesh::mapMesh(const polyMeshMap& map)
{
    forAllIter
    (
        HashPtrTable<globalPolyPatch>,
        patches_,
        iter
    )
    {
        iter()->mapMesh(map);
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

    // Update any patches that have already been added
    polyMesh& mesh = this->db().time().lookupObjectRef<polyMesh>(region);
    forAll(mesh.boundaryMesh(), patchi)
    {
        if (isGlobal(mesh.boundaryMesh()[patchi]))
        {
            patches_
            [
                mesh.boundaryMesh()[patchi].name()
            ]->setDisplacementField(name);
        }
    }
}


void Foam::globalPolyBoundaryMesh::setInverseDisplacement
(
    const word& region,
    const bool inv
)
{
    if (inverseDisplacement_.found(region))
    {
        inverseDisplacement_[region] = inv;
    }
    else
    {
        inverseDisplacement_.insert(region, inv);
    }

    // Update any patches that have already been added
    polyMesh& mesh = this->db().time().lookupObjectRef<polyMesh>(region);
    forAll(mesh.boundaryMesh(), patchi)
    {
        if (isGlobal(mesh.boundaryMesh()[patchi]))
        {
            patches_
            [
                mesh.boundaryMesh()[patchi].name()
            ]->setInverseDisplacement(inv);
        }
    }
}


void Foam::globalPolyBoundaryMesh::clearOut()
{
    forAllIter
    (
        HashPtrTable<globalPolyPatch>,
        patches_,
        iter
    )
    {
        iter()->clearOut();
    }
}


bool Foam::globalPolyBoundaryMesh::write() const
{
    bool good = true;
    forAllIter
    (
        HashPtrTable<globalPolyPatch>,
        patches_,
        iter
    )
    {
        good = good && iter()->write();
    }
    return good;
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
        if (interfaceDicts_.found(mesh().name()))
        {
            const dictionary& mDict = interfaceDicts_[mesh().name()];
            if (mDict.found(pp.name()))
            {
                dict = mDict.subDict(pp.name());
            }
        }
        if (displacementFields_.found(mesh().name()))
        {
            if (!dict.found("displacementField"))
            {
                dict.add
                (
                    "displacementField",
                    displacementFields_[mesh().name()]
                );
            }
        }
        patches_.insert
        (
            pp.name(),
            globalPolyPatch::New(dict, pp).ptr()
        );

        if (inverseDisplacement_.found(mesh().name()))
        {
            patches_[pp.name()]->setInverseDisplacement
            (
                inverseDisplacement_[mesh().name()]
            );
        }
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
    if (!interfaceDicts_.size())
    {
        typeIOobject<IOdictionary> regionPropertiesIO
        (
            IOobject
            (
                "regionProperties",
                mesh().time().constant(),
                mesh().time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
        if (mesh().foundObject<IOdictionary>(regionPropertiesIO.name()))
        {
            const IOdictionary& regionProperties =
                mesh().lookupObject<IOdictionary>(regionPropertiesIO.name());

            if (regionProperties.found("interfaces"))
            {
                interfaceDicts_ =
                    HashTable<dictionary>(regionProperties.lookup("interfaces"));
                readFromRP_ = true;
            }
        }
        else if (regionPropertiesIO.headerOk())
        {
            IOdictionary regionProperties(regionPropertiesIO);
            if (regionProperties.found("interfaces"))
            {
                interfaceDicts_ =
                    HashTable<dictionary>(regionProperties.lookup("interfaces"));
                readFromRP_ = true;
            }
        }
        else if (isA<fvMesh>(mesh()))
        {
            const fvSchemes& schemes =
                dynamicCast<const fvMesh>(mesh()).schemes();
            const entry& e =
                schemes.dict().subDict("interpolationSchemes").lookupEntry
                (
                    pp.name(),
                    false,
                    false
                );
            interfaceDicts_(mesh().name()).set(pp.name(), e.dict());
        }
    }

    bool missingInterpolation = false;
    if (!interfaceDicts_.size())
    {
        if (readFromRP_)
        {
            FatalErrorInFunction
                << "The interfaces is empty in regionProperties. "
                << "This is the default, but a list of "
                << "interfaces is necessary when using coupled patches."
                << "Please specify the interfaces and their mapping methods."
                << "i.e. " << nl
                << "interfaces" << nl
                << "(" << nl
                << "    " << mesh().name()<< nl
                << "    {" << nl
                << "        " << pp.name() << nl
                << "        {" << nl
                << "            ..." << nl
                << "        }" << nl
                << "    }" << nl
                << ");" << endl
                << abort(FatalError);
        }
        else
        {
            missingInterpolation = true;
        }
    }
    if (!interfaceDicts_.found(mesh().name()))
    {
        if (readFromRP_)
        {
            FatalErrorInFunction
                << mesh().name() << " was not found in the list of "
                << "interfaces but a coupled patch was requested for the "
                << "region." << nl
                << "Please specify the region and interface mapping methods"
                << endl
                << abort(FatalError);
        }
        else
        {
            missingInterpolation = true;
        }
    }
    else if (!interfaceDicts_[mesh().name()].isDict(pp.name()))
    {
        if (readFromRP_)
        {
            FatalErrorInFunction
                << pp.name() << " was not found in the list of interfaces "
                << "for region " << mesh().name() << " "
                << "but a coupled patch was requested. Please specify the "
                << "mapping method for the patch" << endl
                << abort(FatalError);
        }
        else
        {
            missingInterpolation = true;
        }
    }

    if (missingInterpolation)
    {
        FatalErrorInFunction
            << "No mapping was provided in fvSchemes/interpolationSchemes. "
            << "Mapping is necessary when using coupled patches."
            << "Please specify the coupled patch mapping methods."
            << "i.e. " << nl
            << "interpolationSchemes" << nl
            << "{" << nl
            << "    " << pp.name() << nl
            << "    {" << nl
            << "        ..." << nl
            << "    }" << nl
            << "}" << endl
            << abort(FatalError);
    }

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
                interfaceDicts_[mesh().name()].subDict(pp.name())
            );
        if (displacementFields_.found(mesh().name()))
        {
            if (!dict.found("displacementField"))
            {
                dict.add
                (
                    "displacementField",
                    displacementFields_[mesh().name()]
                );
            }
        }
        patches_.insert
        (
            pp.name(),
            coupledGlobalPolyPatch::New(dict, pp).ptr()
        );
        if (inverseDisplacement_.found(mesh().name()))
        {
            patches_[pp.name()]->setInverseDisplacement
            (
                inverseDisplacement_[mesh().name()]
            );
        }
    }


    return dynamicCast<const coupledGlobalPolyPatch>(*patches_[pp.name()]);
}


const Foam::coupledGlobalPolyPatch&
Foam::globalPolyBoundaryMesh::operator()(const pointPatch& pp) const
{
    return this->operator()(mesh().boundaryMesh()[pp.name()]);
}


// ************************************************************************* //
