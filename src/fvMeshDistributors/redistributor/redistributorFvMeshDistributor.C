/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "redistributorFvMeshDistributor.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "polyDistributionMap.H"

#include "preserveFaceZonesConstraint.H"
#include "singleProcessorFaceSetsConstraint.H"
#include "preservePatchesConstraint.H"
#include "preserveBafflesConstraint.H"

#include "internalPolyPatch.H"
#include "processorPolyPatch.H"

#include "polyMeshHexRefiner.H"
#include "polyMeshPolyRefiner.H"

#include "hexRefRefinementHistoryConstraint.H"
#include "polyRefinementConstraint.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshDistributors
{
    defineTypeNameAndDebug(redistributor, 0);
    addToRunTimeSelectionTable
    (
        fvMeshDistributor,
        redistributor,
        fvMesh
    );
}
}

Foam::dictionary& Foam::fvMeshDistributors::redistributor::constraints
(
    dictionary& decompositionDict
)
{
    if (!decompositionDict.found("constraints"))
    {
        decompositionDict.set("constraints", dictionary());
    }
    return decompositionDict.subDict("constraints");
}


void Foam::fvMeshDistributors::redistributor::preserveFaceZone
(
    const wordRe& zoneName,
    dictionary& decompositionDict
)
{
    using namespace Foam::decompositionConstraints;
    dictionary& constraintsDict = constraints(decompositionDict);
    wordList toc(constraintsDict.toc());
    forAll(toc, i)
    {
        if (!constraintsDict.isDict(toc[i])) continue;

        dictionary& dict(constraintsDict.subDict(toc[i]));
        word type(dict.lookupOrDefault<word>("type", "none"));

        if (type == preserveFaceZonesConstraint::typeName)
        {
            wordReList zones(dict.lookup("zones"));
            if (!wordReListMatcher(zones).match(zoneName))
            {
                zones.append(zoneName);
                dict.set("zones", zones);
            }
            return;
        }
    }

    constraintsDict.add("faceZones", dictionary());
    dictionary& preserveFaceZonesDict = constraintsDict.subDict("faceZones");
    preserveFaceZonesDict.set
    (
        "type",
        preserveFaceZonesConstraint::typeName
    );
    preserveFaceZonesDict.set("zones", wordReList(1, zoneName));
}


void Foam::fvMeshDistributors::redistributor::singleProcessorFaceSet
(
    const word& setName,
    const label proc,
    dictionary& decompositionDict
)
{
    using namespace Foam::decompositionConstraints;
    dictionary& constraintsDict = constraints(decompositionDict);
    wordList toc(constraintsDict.toc());
    forAll(toc, i)
    {
        if (!constraintsDict.isDict(toc[i])) continue;

        dictionary& dict(constraintsDict.subDict(toc[i]));
        word type(dict.lookupOrDefault<word>("type", "none"));

        if (type == singleProcessorFaceSetsConstraint::typeName)
        {
            List<Tuple2<word, label>> setNameAndProcs
            (
                dict.lookup("singleProcessorFaceSets")
            );
            forAll(setNameAndProcs, i)
            {
                if
                (
                    setName == setNameAndProcs[i].first()
                 && proc == setNameAndProcs[i].second()
                )
                {
                    return;
                }
                else if (setName == setNameAndProcs[i].first())
                {
                    setNameAndProcs[i].second() = proc;
                    dict.set("singleProcessorFaceSets", setNameAndProcs);
                    return;
                }
            }

            setNameAndProcs.append(Tuple2<word, label>(setName, proc));
            dict.set("singleProcessorFaceSets", setNameAndProcs);
            return;
        }
    }

    constraintsDict.add("singleProcessorFaceSets", dictionary());
    dictionary& singleProcessorFaceSetsDict =
        constraintsDict.subDict("singleProcessorFaceSets");
    singleProcessorFaceSetsDict.set
    (
        "type",
        singleProcessorFaceSetsConstraint::typeName
    );
    singleProcessorFaceSetsDict.set
    (
        "singleProcessorFaceSets",
        Tuple2<word, label>(setName, proc)
    );
}


void Foam::fvMeshDistributors::redistributor::preservePatch
(
    const wordRe& patchName,
    dictionary& decompositionDict
)
{
    using namespace Foam::decompositionConstraints;
    dictionary& constraintsDict = constraints(decompositionDict);
    wordList toc(constraintsDict.toc());
    forAll(toc, i)
    {
        if (!constraintsDict.isDict(toc[i])) continue;

        dictionary& dict(constraintsDict.subDict(toc[i]));
        word type(dict.lookupOrDefault<word>("type", "none"));

        if (type == preservePatchesConstraint::typeName)
        {
            wordReList patches(dict.lookup("patches"));
            if (!wordReListMatcher(patches).match(patchName))
            {
                patches.append(patchName);
                dict.set("patches", patches);
            }
            return;
        }
    }

    constraintsDict.add("preservePatches", dictionary());
    dictionary& preservePatchesDict =
        constraintsDict.subDict("preservePatches");
    preservePatchesDict.set
    (
        "type",
        preservePatchesConstraint::typeName
    );
    preservePatchesDict.set("patches", wordReList(1, patchName));
}


void Foam::fvMeshDistributors::redistributor::preserveBaffles
(
    dictionary& decompositionDict
)
{
    using namespace Foam::decompositionConstraints;
    dictionary& constraintsDict = constraints(decompositionDict);
    wordList toc(constraintsDict.toc());
    forAll(toc, i)
    {
        if (!constraintsDict.isDict(toc[i])) continue;

        dictionary& dict(constraintsDict.subDict(toc[i]));
        word type(dict.lookupOrDefault<word>("type", "none"));

        if (type == preserveBafflesConstraint::typeName)
        {
            return;
        }
    }
    constraintsDict.add("preserveBaffles", dictionary());
    dictionary& preserveBafflesDict =
        constraintsDict.subDict("preserveBaffles");
    preserveBafflesDict.set
    (
        "type",
        preserveBafflesConstraint::typeName
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fvMeshDistributors::redistributor::readDict()
{
    if (Pstream::parRun())
    {
        readDict(this->dict());
        return;
    }
}


void Foam::fvMeshDistributors::redistributor::readDict
(
    const dictionary& balanceDict
)
{
    if (!Pstream::parRun())
    {
        balance_ = false;
        return;
    }

    balance_ = balanceDict.lookupOrDefault("balance", true);

    if (!balance_)
    {
        return;
    }

    // Change decomposition method if entry is present
    if (balanceDict.found("method") || balanceDict.found("decomposer"))
    {
        word method =
            balanceDict.lookupBackwardsCompatible
            (
                {"decomposer", "method"}
            );

        bool update = !decomp_.valid();
        if (!update)
        {
            update = decomp_->type() != method;
        }

        if (update)
        {
            decompositionDict_.set("decomposer", method);

            // Remove optional coeffs dictionary since it would override entries and
            // not necessarily be overridden
            if (decompositionDict_.isDict(method + "Coeffs"))
            {
                decompositionDict_.remove(method + "Coeffs");
            }
            decompositionDict_ <<= balanceDict;
            decomp_.clear();
        }

        if (balanceDict.found("contraints"))
        {
            const dictionary& constraintsDict = balanceDict.subDict("constraints");
            forAllConstIter(dictionary, constraintsDict, iter)
            {
                const entry& e = *iter;
                if (e.isDict())
                {
                    addConstraint(e.keyword(), e.dict());
                    decomp_.clear();
                }
            }
        }
    }

    balanceDict.readIfPresent("force", force_);
    balanceDict.readIfPresent("balanceInterval", balanceInterval_);
    balanceDict.readIfPresent("redistributionInterval", balanceInterval_);

    balanceDict.readIfPresent("beginBalance", beginBalance_);
    balanceDict.readIfPresent("endBalance", endBalance_);

    balanceDict.readIfPresent("maxImbalance", maxImbalance_);
}


bool Foam::fvMeshDistributors::redistributor::constraintFound
(
    const word& type
) const
{
    const dictionary& constraintDict = decompositionDict_.subDict("constraints");
    forAllConstIter(dictionary, constraintDict, iter)
    {
        if (iter().isDict())
        {
            if (iter().dict().lookup<word>("type") == type)
            {
                return true;
            }
        }
    }
    return false;
}


void Foam::fvMeshDistributors::redistributor::addConstraint
(
    const word& name,
    const dictionary& dict
)
{
    dictionary& constraintDict = decompositionDict_.subDict("constraints");
    constraintDict.set(name, dict);
}


Foam::decompositionMethod& Foam::fvMeshDistributors::redistributor::decomp()
{
    if (!decomp_.valid())
    {
        decomp_ = decompositionMethod::NewDistributor(decompositionDict_);
    }
    return decomp_();
}


void Foam::fvMeshDistributors::redistributor::distribute
(
    const labelList& distribution,
    const bool dist
)
{
    fvMesh& mesh = this->mesh();

    // Mesh distribution engine
    fvMeshDistribute distributor(mesh);

    // Do actual sending/receiving of mesh
    autoPtr<polyDistributionMap> map
    (
        distributor.distribute(distribution)
    );

    if (dist)
    {
        // Distribute the mesh data
        mesh.distribute(map);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshDistributors::redistributor::redistributor
(
    fvMesh& mesh,
    const bool read
)
:
    fvMeshDistributor(mesh),
    decompositionDict_(decompositionMethod::decomposeParDict(mesh.time())),
    decomp_(nullptr),
    balance_(true),
    force_(false),
    balanceInterval_(10),
    beginBalance_(0),
    endBalance_(great),
    maxImbalance_(0.1),
    timeIndex_(-1),
    iter_(0)
{
    {
        typeIOobject<IOdictionary> dictHeader
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh.time().constant(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        );

        if (dictHeader.headerOk())
        {
            IOdictionary dict(dictHeader);

            if (dict.isDict("distributor"))
            {
                const dictionary& distributorDict = dict.subDict("distributor");

                if (distributorDict.isDict("decomposition"))
                {
                    decompositionDict_ = distributorDict.subDict("decomposition");
                }
            }
        }
    }

    if
    (
        mesh.foundObject<polyMeshHexRefiner>(polyMeshHexRefiner::typeName)
     && !constraintFound(hexRefRefinementHistoryConstraint::typeName)
    )
    {
        // Added refinement history decomposition constraint to keep all
        // cells with the same parent together
        dictionary refinementHistoryDict("refinementHistory");
        refinementHistoryDict.add
        (
            "type",
            hexRefRefinementHistoryConstraint::typeName
        );
        addConstraint("refinementHistory", refinementHistoryDict);
    }
    else if
    (
        mesh.foundObject<polyMeshPolyRefiner>(polyMeshPolyRefiner::typeName)
     && !constraintFound(polyRefinementConstraint::typeName)
    )
    {
        // Added refinement history decomposition constraint to keep all
        // cells with the same parent together
        dictionary refinementHistoryDict("refinementHistory");
        refinementHistoryDict.add
        (
            "type",
            polyRefinementConstraint::typeName
        );
        addConstraint("refinementHistory", refinementHistoryDict);
    }

    if (read)
    {
        readDict();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshDistributors::redistributor::~redistributor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshDistributors::redistributor::update()
{
    const fvMesh& mesh = this->mesh();

    // Only balance if first time step or if some topological changes
    // have happened
    if
    (
        balance_
     && (
            force_
         || (
                Pstream::nProcs() > 1
             && timeIndex_ != mesh.time().timeIndex()
             && (
                    timeIndex_ == 0
                 || (
                        mesh.topoChanging()
                     && ++iter_ % balanceInterval_ == 0
                    )
                )
            )
        )
    )
    {
        timeIndex_ = mesh.time().timeIndex();

        const scalar idealNCells =
            mesh.globalData().nTotalCells()/Pstream::nProcs();

        const scalar imbalance = mag(1 - mesh.nCells()/idealNCells);

        const scalar globalImbalance = returnReduce(imbalance, maxOp<scalar>());

        Info<<"Maximum imbalance = " << 100*globalImbalance << " %" << endl;

        if (debug)
        {
            Pout<< "Current local imbalance = "
                << 100.0*imbalance << "%, "
                << "nCells = " << mesh.nCells()
                << endl;
        }


        if (globalImbalance < maxImbalance_)
        {
            DebugInfo<< "Current imbalance under limit" << endl;
            return false;
        }

        // Create new decomposition distribution
        const labelList distribution
        (
            decomp().decompose(mesh, scalarField())
        );

        // Check if distribution will improve anything
        labelList procLoadNew(Pstream::nProcs(), 0);
        forAll(distribution, celli)
        {
            procLoadNew[distribution[celli]]++;
        }
        reduce(procLoadNew, sumOp<labelList>());
        if (min(procLoadNew) == 0)
        {
            DebugInfo
                << "New distribtion results in a load of 0. Skipping" << endl;
            return false;
        }
        scalar averageLoadNew
        (
            scalar(sum(procLoadNew))/scalar(Pstream::nProcs())
        );
        scalar maxDevNew(max(mag(procLoadNew - averageLoadNew))/averageLoadNew);

        if (maxDevNew > maxImbalance_*0.99)
        {
            Info
                << "    Not balancing because the new distribution does" << nl
                << "    not improve the load. Skipping" << nl
                << "    old imbalance: " << maxImbalance_ << nl
                << "    new imbalance: " << maxDevNew << nl
                << endl;
            return false;
        }

        Info<< "Redistributing mesh with new imbalance = "
            << 100.0*maxDevNew << endl;
        distribute(distribution, true);

        return true;
    }

    return false;
}


bool Foam::fvMeshDistributors::redistributor::forceUpdate(const bool dist)
{
    const fvMesh& mesh = this->mesh();

    const scalar idealNCells =
        mesh.globalData().nTotalCells()/Pstream::nProcs();

    const scalar imbalance = mag(1 - mesh.nCells()/idealNCells);

    const scalar globalImbalance = returnReduce(imbalance, maxOp<scalar>());

    Info<<"Maximum imbalance = " << 100*globalImbalance << " %" << endl;

    if (debug)
    {
        Pout<< "Current local imbalance = "
            << 100.0*imbalance << "%, "
            << "nCells = " << mesh.nCells()
            << endl;
    }


    if (globalImbalance < maxImbalance_)
    {
        DebugInfo<< "Current imbalance under limit" << endl;
        return false;
    }

    // Create new decomposition distribution
    const labelList distribution
    (
        decomp().decompose(mesh, scalarField())
    );

    // Check if distribution will improve anything
    labelList procLoadNew(Pstream::nProcs(), 0);
    forAll(distribution, celli)
    {
        procLoadNew[distribution[celli]]++;
    }
    reduce(procLoadNew, sumOp<labelList>());
    if (min(procLoadNew) == 0)
    {
        DebugInfo
            << "New distribtion results in a load of 0. Skipping" << endl;
        return false;
    }
    scalar averageLoadNew
    (
        scalar(sum(procLoadNew))/scalar(Pstream::nProcs())
    );
    scalar maxDevNew(max(mag(procLoadNew - averageLoadNew))/averageLoadNew);

    if (maxDevNew > maxImbalance_*0.99)
    {
        Info
            << "    Not balancing because the new distribution does" << nl
            << "    not improve the load. Skipping" << nl
            << "    old imbalance: " << maxImbalance_ << nl
            << "    new imbalance: " << maxDevNew << nl
            << endl;
        return false;
    }

    Info<< "Redistributing mesh with new imbalance = "
        << 100.0*maxDevNew << endl;
    distribute(distribution, dist);

    return true;
}


void Foam::fvMeshDistributors::redistributor::topoChange(const polyTopoChangeMap&)
{}


void Foam::fvMeshDistributors::redistributor::mapMesh(const polyMeshMap&)
{}


void Foam::fvMeshDistributors::redistributor::distribute
(
    const polyDistributionMap&
)
{}


bool Foam::fvMeshDistributors::redistributor::write(const bool write) const
{
    return true;
}


// ************************************************************************* //
