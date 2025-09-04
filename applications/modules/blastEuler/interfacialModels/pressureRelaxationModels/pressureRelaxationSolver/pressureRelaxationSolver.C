/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022-2025
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
}


Foam::autoPtr<Foam::pressureRelaxationSolver>
Foam::pressureRelaxationSolver::New
(
    phaseSystem& fluid,
    interfacialPressureModelTable& interfacialPressureModels,
    pressureRelaxationModelTable& pressureRelaxationModels
)
{
    label nFluids = 0;
    forAll(fluid.phases(), phasei)
    {
        if (!fluid.phases()[phasei].slavePressure())
        {
            nFluids++;
        }
    }
    if (nFluids < 2)
    {
        return autoPtr<pressureRelaxationSolver>
        (
            new pressureRelaxationSolver(fluid, false)
        );
    }
    const dictionary& solverDict = fluid.subDict("pressureRelaxationSolver");
    const word relaxationType = solverDict.lookup<word>("type");
    Info<< "Selecting " << typeName <<  ": " << relaxationType << endl;

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

    return cstrIter()
    (
        solverDict.optionalSubDict(relaxationType + "Coeffs"),
        fluid,
        interfacialPressureModels,
        pressureRelaxationModels
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureRelaxationSolver::pressureRelaxationSolver
(
    phaseSystem& fluid,
    interfacialPressureModelTable& interfacialPressureModels,
    pressureRelaxationModelTable& pressureRelaxationModels
)
:
    pressureRelaxationSolver(fluid, interfacialPressureModels)
{
    forAll(includedPhases_, i)
    {
        const phaseModel& phaseI = fluid.phases()[includedPhases_[i]];
        for (label j = i+1; j < includedPhases_.size(); j++)
        {
            const phaseModel& phaseJ = fluid.phases()[includedPhases_[j]];
            phasePairKey key(phaseI.name(), phaseJ.name());

            if (!pressureRelaxationModels.found(key))
            {
                FatalErrorInFunction
                    << "Did not find pressureRelaxationModel for "
                    << key << endl
                    << abort(FatalError);
            }
            pressureRelaxationModels_.append
            (
                &pressureRelaxationModels[key]()
            );
        }
    }
}


Foam::pressureRelaxationSolver::pressureRelaxationSolver
(
    phaseSystem& fluid,
    interfacialPressureModelTable& interfacialPressureModels
)
:
    pressureRelaxationSolver(fluid, true)
{
    //- Add unorded phase pairs with vaild pressureRelaxation models
    forAll(includedPhases_, i)
    {
        const phaseModel& phaseI = fluid.phases()[includedPhases_[i]];
        for (label j = i+1; j < includedPhases_.size(); j++)
        {
            const phaseModel& phaseJ = fluid.phases()[includedPhases_[j]];
            phasePairKey key(phaseI.name(), phaseJ.name());

            if (!interfacialPressureModels.found(key))
            {
                FatalErrorInFunction
                    << "Did not find interfacialPressureModel for "
                    << key << endl
                    << abort(FatalError);
            }
            interfacialPressureModels_.append
            (
                &interfacialPressureModels[key]()
            );
        }
    }
}


Foam::pressureRelaxationSolver::pressureRelaxationSolver
(
    phaseSystem& fluid,
    const bool needPhases
)
:
    solvePressureRelaxation_(false),
    fluid_(fluid),
    includedPhases_(0),
    phaseModels_(0),
    fixedPhaseModels_(0),
    phaseIndicies_(0),
    thermos_(0),
    interfacialPressureModels_(0),
    pressureRelaxationModels_(0),
    nEqns_(0)
{
    if (!needPhases)
    {
        return;
    }

    forAll(fluid.phases(), phasei)
    {
        if (!fluid.phases()[phasei].slavePressure())
        {
            includedPhases_.append(phasei);
        }
        else
        {
            fixedPhaseModels_.append(&fluid.phases()[phasei]);
        }
    }
    if (includedPhases_.size() <= 1)
    {
        return;
    }

    solvePressureRelaxation_ = true;
    phaseModels_.setSize(includedPhases_.size());
    thermos_.setSize(includedPhases_.size());
    residualAlphas_.setSize(includedPhases_.size());

    hashedWordList includedPhases;
    forAll(includedPhases_, phasei)
    {
        const phaseModel& phase = fluid.phases()[includedPhases_[phasei]];
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
            &fluid.phases()[includedPhases_[phasei]]
        );
        thermos_.set
        (
            phasei,
            &fluid.mesh().lookupObjectRef<fluidBlastThermo>
            (
                IOobject::groupName(physicalProperties::typeName, phase.group())
            )
        );
        residualAlphas_[phasei] = phaseModels_[phasei].residualAlpha().value();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pressureRelaxationSolver::~pressureRelaxationSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
