/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.
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

#include "pressureRelaxationSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pressureRelaxationSolver, 0);
    defineRunTimeSelectionTable(pressureRelaxationSolver, dictionary);
    addNamedToRunTimeSelectionTable
    (
        pressureRelaxationSolver,
        pressureRelaxationSolver,
        dictionary,
        none
    );
}


Foam::autoPtr<Foam::pressureRelaxationSolver>
Foam::pressureRelaxationSolver::New
(
    phaseSystem& fluid,
    interfacialPressureModelTable& interfacialPressureModels,
    pressureRelaxationModelTable& pressureRelaxationModels
)
{
    word relaxationType("none");

    label nFluids = 0;
    forAll(fluid.phases(), phasei)
    {
        if (!fluid.phases()[phasei].slavePressure())
        {
            nFluids++;
        }
    }
    if (nFluids > 1)
    {
        relaxationType = fluid.lookup<word>("pressureRelaxationSolver");
        Info<< "Selecting pressureRelaxationSolver: "
            << relaxationType << endl;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(relaxationType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown pressureRelaxationSolver type "
            << relaxationType << endl << endl
            << "Valid pressureRelaxationSolver types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(fluid, interfacialPressureModels, pressureRelaxationModels);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureRelaxationSolver::pressureRelaxationSolver
(
    phaseSystem& fluid,
    interfacialPressureModelTable& interfacialPressureModels,
    pressureRelaxationModelTable& pressureRelaxationModels,
    const bool needIPModels,
    const bool needPRModels
)
:
    solvePressureRelaxation_(false),
    fluid_(fluid),
    phaseModels_(),
    phaseIndicies_(),
    thermos_(),
    interfacialPressureModels_(interfacialPressureModels.size()),
    pressureRelaxationModels_(pressureRelaxationModels.size()),
    nEqns_(0)
{
    labelList phases;
    forAll(fluid.phases(), phasei)
    {
        if (!fluid.phases()[phasei].slavePressure())
        {
            phases.append(phasei);
        }
    }
    if (phases.size() <= 1)
    {
        return;
    }

    solvePressureRelaxation_ = true;
    phaseModels_.setSize(phases.size());
    thermos_.setSize(phases.size());

    hashedWordList includedPhases;
    forAll(phases, phasei)
    {
        const phaseModel& phase = fluid.phases()[phases[phasei]];
        if (phase.slavePressure())
        {
            FatalErrorInFunction
                << "Trying to relax " << phase.name()
                << " but it has a slave pressure" << endl
                << abort(FatalError);
        }
        includedPhases.append(phase.name());
        phaseIndicies_.insert(phase.index(), phasei);
        phaseModels_.set
        (
            phasei,
            &fluid.phases()[phases[phasei]]
        );
        thermos_.set
        (
            phasei,
            &fluid.mesh().lookupObjectRef<fluidBlastThermo>
            (
                IOobject::groupName(basicThermo::dictName, phase.group())
            )
        );
    }

    //- Add unorded phase pairs with vaild pressureRelaxation models
    if (needIPModels)
    {
        label pairi = 0;
        forAll(phases, i)
        {
            const phaseModel& phaseI = fluid.phases()[phases[i]];
            for (label j = i+1; j < phases.size(); j++)
            {
                const phaseModel& phaseJ = fluid.phases()[phases[j]];
                phasePairKey key(phaseI.name(), phaseJ.name());

                if (!interfacialPressureModels.found(key))
                {
                    FatalErrorInFunction
                        << "Did not find interfacialPressureModel for "
                        << key << endl
                        << abort(FatalError);
                }
                interfacialPressureModels_.set
                (
                    pairi++,
                    &interfacialPressureModels[key]()
                );
            }
        }
        interfacialPressureModels_.resize(pairi);
    }

    if (needPRModels)
    {
        label pairi = 0;
        forAll(phases, i)
        {
            const phaseModel& phaseI = fluid.phases()[phases[i]];
            for (label j = i+1; j < phases.size(); j++)
            {
                const phaseModel& phaseJ = fluid.phases()[phases[j]];
                phasePairKey key(phaseI.name(), phaseJ.name());

                if (!pressureRelaxationModels.found(key))
                {
                    FatalErrorInFunction
                        << "Did not find pressureRelaxationModel for "
                        << key << endl
                        << abort(FatalError);
                }
                pressureRelaxationModels_.set
                (
                    pairi++,
                    &pressureRelaxationModels[key]()
                );
            }
        }
    }
}


Foam::pressureRelaxationSolver::pressureRelaxationSolver
(
    phaseSystem& fluid,
    interfacialPressureModelTable& interfacialPressureModels,
    pressureRelaxationModelTable& pressureRelaxationModels
)
:
    solvePressureRelaxation_(false),
    fluid_(fluid),
    phaseModels_(0),
    phaseIndicies_(0),
    thermos_(0),
    interfacialPressureModels_(0),
    pressureRelaxationModels_(0),
    nEqns_(0)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pressureRelaxationSolver::~pressureRelaxationSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
