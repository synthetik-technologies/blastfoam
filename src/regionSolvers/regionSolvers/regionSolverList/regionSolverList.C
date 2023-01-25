/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

#include "regionSolverList.H"
#include "regionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionSolverList, 0);
}


// * * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolverList::regionSolverList
(
    const dictionary& regionProperties,
    const Time& runTime,
    PtrList<dynamicFvMesh>& regionMeshes,
    const List<Pair<word>>& regionTypes
)
:
    PtrList<regionSolver>(regionMeshes.size()),
    runTime_(runTime),
    regionProperties_(regionProperties),
    changed_(this->size(), true),
    predictSolids_(solutionControls().lookupOrDefault("predictSolids", true)),
    fixedMapping_(regionProperties_.lookupOrDefault<bool>("fixedMapping", false)),
    iterNo_(0)
{
    // Clearing of global patches is handled internally
    globalPolyBoundaryMesh::clearOnMovement = false;
    forAll(regionMeshes, regioni)
    {
        this->set
        (
            regioni,
            regionSolver::New
            (
                regionTypes[regioni].second(),
                regionMeshes[regioni],
                *this
            ).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolverList::~regionSolverList()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::regionSolverList::converged() const
{
    bool allUnknown = true;
    forAll(*this, regioni)
    {
        switch (operator[](regioni).converged())
        {
            case NOT_CONVERGED:
                return false;
            case CONVERGED:
                allUnknown = false;
                break;
            default:
                break;
        }
    }
    return !allUnknown;
}


Foam::Convergence Foam::regionSolverList::convergence() const
{
    bool allUnknown = true;
    forAll(*this, regioni)
    {
        switch (operator[](regioni).converged())
        {
            case NOT_CONVERGED:
                return NOT_CONVERGED;
            case CONVERGED:
                allUnknown = false;
                break;
            default:
                break;
        }
    }
    return allUnknown ? UNKNOWN_CONVERGENCE : CONVERGED;
}

void Foam::regionSolverList::initialiseDisplacement()
{
    // Make sure all regions are setup
    forAll(*this, regioni)
    {
        operator[](regioni).initialiseFields();
    }

    if (debug)
    {
        const_cast<Time&>(runTime_).setTime(runTime_.value()+runTime_.deltaTValue(), runTime_.timeIndex());
    }

    update(true);

    // Initialize meshes
    iterNo_ = 0;
    Info<< endl;
    IOobject::writeDivider(Info)
        << "Initial unrelaxed displacement iteration" << nl << endl;
    forAll(*this, regioni)
    {
        Info<< operator[](regioni).name() << ": " << endl;
        operator[](regioni).initialiseMesh(FIRST_ITER);
        Info<< endl;
    }

    // Set displacement field names
    initialise();

    label nInitialCorrectors =
        solutionControls().lookup<label>("nInitialCorrectors");

    do
    {
        if (debug)
        {
            const_cast<Time&>(runTime_).writeNow();
            const_cast<Time&>(runTime_).setTime(runTime_.value()+runTime_.deltaTValue(), runTime_.timeIndex());
        }

        Info<< endl;
        IOobject::writeDivider(Info)
            << "Initial correction iteration: " << iterNo_ << nl << endl;

        // Force updating of mapping
        update(true);

        // Initialize meshes
        forAll(*this, regioni)
        {
            operator[](regioni).moveMesh(MID_ITER);
        }
    } while (iterNo_++ < nInitialCorrectors && !converged());

    if (debug)
    {
        const_cast<Time&>(runTime_).writeNow();
        const_cast<Time&>(runTime_).setTime(runTime_.value()+runTime_.deltaTValue(), runTime_.timeIndex());

    }

    if (convergence() == CONVERGED)
    {
        Info<< "Converged initial displacement in " << iterNo_
            << " iterations" << nl << endl;
    }
    else if (convergence() == NOT_CONVERGED)
    {
        Info<< "*** Initial displacement did not converge" << nl << endl;
    }

    // Final update with no relaxation
    Info<< "Final unrelaxed iteration" << nl << endl;
    forAll(*this, regioni)
    {
        operator[](regioni).initialiseMesh(FINAL_ITER);
    }

    if (debug)
    {
        const_cast<Time&>(runTime_).writeNow();
    }
}


void Foam::regionSolverList::update(const bool force)
{
    forAll(*this, regioni)
    {
        if (changed_[regioni] || force)
        {
            operator[](regioni).update();
        }
    }
    changed_ = false;
}


bool Foam::regionSolverList::changeMesh()
{
    forAll(*this, regioni)
    {
        changed_[regioni] = operator[](regioni).changeMesh();
    }
    bool changed = anyChanged();

    // Mapping ALWAYS needs to be updated if there have been changes
    // Update immediately
    update();

    return changed;
}


void Foam::regionSolverList::initialiseMesh(const IterType iter)
{
    forAll(*this, regioni)
    {
        operator[](regioni).initialiseMesh(iter);
    }
}


void Foam::regionSolverList::initialiseFields()
{
    forAll(*this, regioni)
    {
        operator[](regioni).initialiseFields();
    }
}


void Foam::regionSolverList::initialise()
{
    forAll(*this, regioni)
    {
        if (operator[](regioni).isSolid())
        {
            operator[](regioni).initialise();
        }
    }
    forAll(*this, regioni)
    {
        if (!operator[](regioni).isSolid())
        {
            operator[](regioni).initialise();
        }
    }
    update(fixedMapping_);
}


bool Foam::regionSolverList::moveMesh(const IterType iter)
{
    forAll(*this, regioni)
    {
        changed_[regioni] =
            changed_[regioni] || operator[](regioni).moveMesh(iter);
    }
    return anyChanged();
}


void Foam::regionSolverList::solve()
{
    label nOuterCorrectors =
        solutionControls().lookup<label>("nOuterCorrectors");
    iterNo_ = 0;
    bool finished = false;
    bool cleanup = false;
    do
    {
        Info<< endl;
        IOobject::writeDivider(Info)
            << "Outer iteration: " << iterNo_ << nl
            << "Time = " << runTime_.timeName() << nl << endl;

        IterType iter =
            iterNo_ == 0 ? FIRST_ITER
          : (
                (cleanup || iterNo_ == nOuterCorrectors-1)
              ? FINAL_ITER
              : MID_ITER
            );
        if (cleanup)
        {
            finished = true;
        }

        if (predictSolids_)
        {
            forAll(*this, regioni)
            {
                if (operator[](regioni).isSolid())
                {
                    Info<< "Solving region "
                        << operator[](regioni).mesh().name() << endl;
                    operator[](regioni).solve();
                }
            }
            forAll(*this, regioni)
            {
                if (!operator[](regioni).isSolid())
                {
                    Info<< "Solving region "
                        << operator[](regioni).mesh().name() << endl;
                    operator[](regioni).moveMesh(iter);
                    operator[](regioni).solve();
                }
            }
        }
        else
        {
            forAll(*this, regioni)
            {
                if (!operator[](regioni).isSolid())
                {
                    Info<< "Solving region "
                        << operator[](regioni).mesh().name() << endl;
                    operator[](regioni).moveMesh(iter);
                    operator[](regioni).solve();
                }
            }
            forAll(*this, regioni)
            {
                if (operator[](regioni).isSolid())
                {
                    Info<< "Solving region "
                        << operator[](regioni).mesh().name() << endl;
                    operator[](regioni).solve();
                }
            }
        }
        iterNo_++;

        if (converged())
        {
            cleanup = true;
        }

    } while (!finished && iterNo_ < nOuterCorrectors);

    if (convergence() == CONVERGED)
    {
        Info<< "All regions converged in " << iterNo_
            << " iterations" << nl << endl;
    }
    else if (convergence() == NOT_CONVERGED)
    {
        Info<< "*** Regions did not converge" << nl << endl;
    }

    update();
    clear();
}


void Foam::regionSolverList::clear()
{
    forAll(*this, regioni)
    {
        operator[](regioni).clear();
    }
}


Foam::scalar Foam::regionSolverList::CoNum() const
{
    scalar co = 0.0;
    forAll(*this, regioni)
    {
        co = max(co, operator[](regioni).CoNum());
    }
    return co;
}


Foam::scalar Foam::regionSolverList::maxCo() const
{
    scalar co = great;
    forAll(*this, regioni)
    {
        co = min(co, operator[](regioni).maxCo());
    }
    return co;
}


Foam::scalar Foam::regionSolverList::newDeltaT() const
{
    scalar deltaT = great;
    forAll(*this, regioni)
    {
        deltaT = min(deltaT, operator[](regioni).newDeltaT());
    }
    Info<< "deltaT = " <<  deltaT << endl;
    return deltaT;
}

// ************************************************************************* //
