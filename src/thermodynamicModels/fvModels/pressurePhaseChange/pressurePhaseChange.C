/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     | Website:  https://openfoam.org
   \\  /    A nd           | Copyright (C) 2019-2026
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

#include "pressurePhaseChange.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"
#include "multiphaseFluidBlastThermo.H"
#include "multicomponentBlastThermo.H"
#include "twoPhaseFluidBlastThermo.H"
#include "physicalProperties.H"
#include "Time.H"
#include "fvMatrices.H"
#include "fvTimeIntegrator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pressurePhaseChange, 0);
addToRunTimeSelectionTable(fvModel, pressurePhaseChange, dictionary);
addNamedToRunTimeSelectionTable(fvModel, pressurePhaseChange, dictionary, pressurePhaseTransfer);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pressurePhaseChange::pressurePhaseChange
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    sourcePhaseName_(),
    targetPhases_(),
    pLower_(0.0),
    pUpper_(1e7),
    pRange_(1e6),
    tau_(1e-4),
    isMoleFraction_(false),
    fractionsConverted_(false),
    massRate_(0.0),
    eOffset_("eOffset", dimEnergy/dimMass, 0.0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> pressurePhaseChange::gamma() const
{
    const volScalarField& p = mesh().lookupObject<volScalarField>("p");

    dimensionedScalar pUpper("pUpper", dimPressure, pUpper_);
    dimensionedScalar pRange("pRange", dimPressure, pRange_);

    if (pLower_ > 0.0)
    {
        dimensionedScalar pLower("pLower", dimPressure, pLower_);
        return 
            0.25 * (scalar(1.0) + tanh((pUpper - p)/pRange))
                 * (scalar(1.0) + tanh((p - pLower)/pRange))
                 * pos(pUpper + scalar(2.0)*pRange - p)
                 * pos(p - (pLower - scalar(2.0)*pRange));
    }
    else
    {
        return 
            0.5 * (scalar(1.0) + tanh((pUpper - p)/pRange))
                * pos(pUpper + scalar(2.0)*pRange - p);
    }
}


void pressurePhaseChange::lookupPhaseThermo
(
    const word& phaseName,
    const phaseBlastThermo*& mixturePhasePtr,
    const blastThermo*& phaseThermoPtr
) const
{
    mixturePhasePtr = nullptr;
    phaseThermoPtr = nullptr;

    // Per-phase thermo (blastEulerFoam)
    const word phaseDictName
    (
        IOobject::groupName(physicalProperties::typeName, phaseName)
    );
    if (mesh().foundObject<fluidBlastThermo>(phaseDictName))
    {
        phaseThermoPtr =
            &mesh().lookupObject<fluidBlastThermo>(phaseDictName);
        return;
    }

    // Phase of the mixture thermo (blastFoam two-/multiphase systems)
    if (!mesh().foundObject<fluidBlastThermo>(physicalProperties::typeName))
    {
        return;
    }

    const fluidBlastThermo& mixture =
        mesh().lookupObject<fluidBlastThermo>(physicalProperties::typeName);

    if (isA<multiphaseFluidBlastThermo>(mixture))
    {
        const multiphaseFluidBlastThermo& thermo =
            refCast<const multiphaseFluidBlastThermo>(mixture);
        forAll(thermo.phaseNames(), phasei)
        {
            if (thermo.thermo(phasei).phaseName() == phaseName)
            {
                mixturePhasePtr = &thermo.thermo(phasei);
                return;
            }
        }
    }
    else if (isA<twoPhaseFluidBlastThermo>(mixture))
    {
        const twoPhaseFluidBlastThermo& thermo =
            refCast<const twoPhaseFluidBlastThermo>(mixture);
        for (label phasei = 0; phasei < 2; phasei++)
        {
            if (thermo.thermo(phasei).phaseName() == phaseName)
            {
                mixturePhasePtr = &thermo.thermo(phasei);
                return;
            }
        }
    }
}


tmp<volScalarField> pressurePhaseChange::energyOffset() const
{
    // Constant energy-datum shift [J/kg] added to the transferred mass, for
    // cases where the source EOS is deliberately kept on a different energy
    // datum from the target phases. A state-dependent offset is not used:
    // it would add or remove real energy rather than relabel the datum.
    // Reconcile the datums through the source EOS E0 instead where possible.
    return volScalarField::New("energyOffset", mesh(), eOffset_);
}


void pressurePhaseChange::convertMoleFractions() const
{
    if (!isMoleFraction_ || fractionsConverted_)
    {
        return;
    }

    for (targetPhase& tp : targetPhases_)
    {
        if (tp.targetFractions.empty())
        {
            continue;
        }

        // Find the multicomponent thermo of the target phase
        const phaseBlastThermo* mixturePhasePtr = nullptr;
        const blastThermo* phaseThermoPtr = nullptr;
        lookupPhaseThermo(tp.name, mixturePhasePtr, phaseThermoPtr);

        const multicomponentBlastThermo* compPtr =
            mixturePhasePtr
          ? dynamic_cast<const multicomponentBlastThermo*>(mixturePhasePtr)
          : dynamic_cast<const multicomponentBlastThermo*>(phaseThermoPtr);

        if (!compPtr)
        {
            FatalErrorInFunction
                << "isMoleFraction is true but no multicomponent thermo "
                << "was found for target phase " << tp.name
                << " to provide species molecular weights"
                << exit(FatalError);
        }

        const speciesTable& species = compPtr->species();

        scalar sumXW = 0;
        forAllConstIter(HashTable<scalar>, tp.targetFractions, iter)
        {
            if (!species.found(iter.key()))
            {
                FatalErrorInFunction
                    << "Specie " << iter.key() << " in targetFractions "
                    << "is not a specie of phase " << tp.name << nl
                    << "Valid species are " << species
                    << exit(FatalError);
            }
            sumXW += iter()*compPtr->WiValue(species[iter.key()]);
        }

        if (sumXW <= 0)
        {
            FatalErrorInFunction
                << "Mole fractions for target phase " << tp.name
                << " sum to zero" << exit(FatalError);
        }

        // Y_i = X_i W_i/sum_j(X_j W_j)
        Info<< "    " << name() << ": target mass fractions for phase "
            << tp.name << " (converted from mole fractions):" << nl;
        forAllIter(HashTable<scalar>, tp.targetFractions, iter)
        {
            iter() *= compPtr->WiValue(species[iter.key()])/sumXW;
            Info<< "        " << iter.key() << " " << iter() << nl;
        }
        Info<< endl;
    }

    fractionsConverted_ = true;
}


wordList pressurePhaseChange::addSupFields() const
{
    if (mesh().foundObject<fvTimeIntegrator>("timeIntegrator"))
    {
        return wordList();
    }

    label nFields = 5; // rho, rhoU, rhoE, rho.source, alphaRho.source
    for (const targetPhase& tp : targetPhases_)
    {
        nFields += 2; // rho.target, alphaRho.target
        nFields += tp.targetFractions.size(); // specie.target
    }

    wordList fields(nFields);
    label i = 0;
    fields[i++] = "rho";
    fields[i++] = "rhoU";
    fields[i++] = "rhoE";
    fields[i++] = "rho." + sourcePhaseName_;
    fields[i++] = "alphaRho." + sourcePhaseName_;

    for (const targetPhase& tp : targetPhases_)
    {
        fields[i++] = "rho." + tp.name;
        fields[i++] = "alphaRho." + tp.name;
        forAllConstIter(HashTable<scalar>, tp.targetFractions, iter)
        {
            fields[i++] = iter.key() + "." + tp.name;
        }
    }

    return fields;
}


bool pressurePhaseChange::addsSupToField(const word& fieldName) const
{
    if (mesh().foundObject<fvTimeIntegrator>("timeIntegrator"))
    {
        return false;
    }

    if (fieldName == "rho" || fieldName == "rhoU" || fieldName == "rhoE" || fieldName == "e" || fieldName == "he")
    {
        return true;
    }

    if (fieldName.find("rho.") == 0 || fieldName.find("alphaRho.") == 0 || fieldName.find("alpha.") == 0)
    {
        return true;
    }

    for (const targetPhase& tp : targetPhases_)
    {
        size_t dotIndex = fieldName.find('.');
        if (dotIndex != string::npos)
        {
            word phaseName = fieldName.substr(dotIndex + 1);
            if (phaseName == tp.name)
            {
                word specieName = fieldName.substr(0, dotIndex);
                if (tp.targetFractions.found(specieName))
                {
                    return true;
                }
            }
        }
    }

    return false;
}


void pressurePhaseChange::preUpdateMesh()
{
    fvModel::preUpdateMesh();

    convertMoleFractions();

    if (mesh().foundObject<fvTimeIntegrator>("timeIntegrator"))
    {
        fvTimeIntegrator& integrator = 
            const_cast<fvTimeIntegrator&>
            (
                mesh().lookupObject<fvTimeIntegrator>("timeIntegrator")
            );

        const volScalarField& alphaSource = 
            mesh().lookupObject<volScalarField>("alpha." + sourcePhaseName_);
        const volScalarField& rhoSource = 
            mesh().lookupObject<volScalarField>("rho." + sourcePhaseName_);

        tmp<volScalarField> tgamma = gamma();
        
        volScalarField mDot = tgamma() * alphaSource * rhoSource / tau_;
        if (mesh().foundObject<volScalarField>("lambda." + sourcePhaseName_))
        {
            const volScalarField& lambda = 
                mesh().lookupObject<volScalarField>("lambda." + sourcePhaseName_);
            mDot *= lambda;
        }
        
        dimensionedScalar residualRho("residualRho", dimDensity, 1e-10);

        // Mass sources [kg/m3/s]
        // rho.Phase
        volScalarField::Internal SrhoSource
        (
            IOobject("SrhoSource", mesh().time().name(), mesh()),
            mesh(),
            dimensionedScalar("SrhoSource", dimDensity/dimTime, 0.0)
        );
        SrhoSource.primitiveFieldRef() = mDot.primitiveField();
        integrator.addSource("rho." + sourcePhaseName_, SrhoSource);

        for (const targetPhase& tp : targetPhases_)
        {
            volScalarField::Internal SrhoTarget
            (
                IOobject("SrhoTarget", mesh().time().name(), mesh()),
                mesh(),
                dimensionedScalar("SrhoTarget", dimDensity/dimTime, 0.0)
            );
            SrhoTarget.primitiveFieldRef() = -mDot.primitiveField() * tp.fraction;
            integrator.addSource("rho." + tp.name, SrhoTarget);

            // alphaRho.Phase
            integrator.addSource("alphaRho." + tp.name, SrhoTarget);

            // Alpha source [1/s]
            if (mesh().foundObject<volScalarField>("alpha." + tp.name))
            {
                volScalarField::Internal SalphaTarget
                (
                    IOobject("SalphaTarget", mesh().time().name(), mesh()),
                    mesh(),
                    dimensionedScalar("SalphaTarget", dimless/dimTime, 0.0)
                );
                SalphaTarget.primitiveFieldRef() = -mDot.primitiveField() * tp.fraction / max(rhoSource.primitiveField(), residualRho.value());
                integrator.addSource("alpha." + tp.name, SalphaTarget);
            }

            // Species sources [kg/m3/s]
            forAllConstIter(HashTable<scalar>, tp.targetFractions, iter)
            {
                volScalarField::Internal SY
                (
                    IOobject("SY", mesh().time().name(), mesh()),
                    mesh(),
                    dimensionedScalar("SY", dimDensity/dimTime, 0.0)
                );
                SY.primitiveFieldRef() = -mDot.primitiveField() * tp.fraction * iter();
                integrator.addSource(iter.key() + "." + tp.name, SY);
            }
        }

        // alphaRho.Phase for source
        integrator.addSource("alphaRho." + sourcePhaseName_, SrhoSource);

        // Energy source for rhoE or per-phase alphaRhoE
        tmp<volScalarField> tEOff = energyOffset();
        const volScalarField& eOff = tEOff();
        if (mesh().foundObject<volScalarField::Internal>("rhoE"))
        {
            volScalarField::Internal SrhoE
            (
                IOobject("SrhoE", mesh().time().name(), mesh()),
                mesh(),
                dimensionedScalar("SrhoE", dimEnergy/dimVolume/dimTime, 0.0)
            );
            SrhoE.primitiveFieldRef() = -mDot.primitiveField() * eOff.primitiveField();
            integrator.addSource("rhoE", SrhoE);
        }
        else
        {
            // Multi-fluid (blastEulerFoam): each phase carries its own
            // momentum and total energy, so the transferred mass takes the
            // source-phase velocity U_s and specific total energy
            // E_s = alphaRhoE_s/alphaRho_s with it. The targets also receive
            // the constant datum shift eOffset. Mass, momentum and
            // (datum-corrected) energy are conserved exactly.
            const word alphaRhoSName("alphaRho." + sourcePhaseName_);
            const word alphaRhoESName("alphaRhoE." + sourcePhaseName_);
            const word USName("U." + sourcePhaseName_);

            if
            (
                !mesh().foundObject<volScalarField>(alphaRhoSName)
             || !mesh().foundObject<volScalarField>(alphaRhoESName)
             || !mesh().foundObject<volVectorField>(USName)
            )
            {
                FatalErrorInFunction
                    << "Cannot find " << alphaRhoSName << ", "
                    << alphaRhoESName << " and " << USName
                    << " for the per-phase energy and momentum transfer"
                    << exit(FatalError);
            }

            const scalarField& alphaRhoS =
                mesh().lookupObject<volScalarField>(alphaRhoSName);
            const scalarField& alphaRhoES =
                mesh().lookupObject<volScalarField>(alphaRhoESName);
            const vectorField& US =
                mesh().lookupObject<volVectorField>(USName);

            const scalarField ES(alphaRhoES/max(alphaRhoS, small));
            const scalarField& mDotI = mDot.primitiveField();
            const scalarField& eOffI = eOff.primitiveField();

            // Source phase: loses mDot*U_s and mDot*E_s
            volScalarField::Internal SalphaRhoESource
            (
                IOobject("SalphaRhoESource", mesh().time().name(), mesh()),
                mesh(),
                dimensionedScalar(dimEnergy/dimVolume/dimTime, 0.0)
            );
            SalphaRhoESource.primitiveFieldRef() = mDotI*ES;
            integrator.addSource(alphaRhoESName, SalphaRhoESource);

            const word alphaRhoUSName("alphaRhoU." + sourcePhaseName_);
            if (mesh().foundObject<volVectorField>(alphaRhoUSName))
            {
                volVectorField::Internal SalphaRhoUSource
                (
                    IOobject("SalphaRhoUSource", mesh().time().name(), mesh()),
                    mesh(),
                    dimensionedVector(dimDensity*dimVelocity/dimTime, Zero)
                );
                SalphaRhoUSource.primitiveFieldRef() = mDotI*US;
                integrator.addSource(alphaRhoUSName, SalphaRhoUSource);
            }

            // Target phases: gain fraction*mDot*U_s and
            // fraction*mDot*(E_s + eOffset)
            for (const targetPhase& tp : targetPhases_)
            {
                const word alphaRhoETName("alphaRhoE." + tp.name);
                if (mesh().foundObject<volScalarField>(alphaRhoETName))
                {
                    volScalarField::Internal SalphaRhoETarget
                    (
                        IOobject
                        (
                            "SalphaRhoETarget",
                            mesh().time().name(),
                            mesh()
                        ),
                        mesh(),
                        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0.0)
                    );
                    SalphaRhoETarget.primitiveFieldRef() =
                        -tp.fraction*mDotI*(ES + eOffI);
                    integrator.addSource(alphaRhoETName, SalphaRhoETarget);
                }

                const word alphaRhoUTName("alphaRhoU." + tp.name);
                if (mesh().foundObject<volVectorField>(alphaRhoUTName))
                {
                    volVectorField::Internal SalphaRhoUTarget
                    (
                        IOobject
                        (
                            "SalphaRhoUTarget",
                            mesh().time().name(),
                            mesh()
                        ),
                        mesh(),
                        dimensionedVector(dimDensity*dimVelocity/dimTime, Zero)
                    );
                    SalphaRhoUTarget.primitiveFieldRef() =
                        -tp.fraction*mDotI*US;
                    integrator.addSource(alphaRhoUTName, SalphaRhoUTarget);
                }
            }
        }

        // Alpha source [1/s]
        if (mesh().foundObject<volScalarField>("alpha." + sourcePhaseName_))
        {
            volScalarField::Internal Salpha
            (
                IOobject("Salpha", mesh().time().name(), mesh()),
                mesh(),
                dimensionedScalar("Salpha", dimless/dimTime, 0.0)
            );
            Salpha.primitiveFieldRef() = (mDot.primitiveField() / max(rhoSource.primitiveField(), residualRho.value()));
            integrator.addSource("alpha." + sourcePhaseName_, Salpha);
        }

        massRate_ = gSum(mDot.primitiveField() * mesh().V());
    }
}


void pressurePhaseChange::correct()
{
    scalar totalMass = massRate_ * mesh().time().deltaTValue();
    Info << "    " << name() << ": total mass transferred = " << totalMass << " kg" << endl;
}


void pressurePhaseChange::addSup
(
    const volScalarField& field,
    fvMatrix<scalar>& eqn
) const
{
    if (mesh().foundObject<fvTimeIntegrator>("timeIntegrator"))
    {
        return;
    }

    convertMoleFractions();

    const word& fieldName = field.name();
    const volScalarField& alphaSource = 
        mesh().lookupObject<volScalarField>("alpha." + sourcePhaseName_);
    const volScalarField& rhoSource = 
        mesh().lookupObject<volScalarField>("rho." + sourcePhaseName_);

    tmp<volScalarField> tgamma = gamma();
    
    volScalarField mDot = tgamma() * alphaSource * rhoSource / tau_;
    if (mesh().foundObject<volScalarField>("lambda." + sourcePhaseName_))
    {
        const volScalarField& lambda = 
            mesh().lookupObject<volScalarField>("lambda." + sourcePhaseName_);
        mDot *= lambda;
    }

    if (fieldName == "rho." + sourcePhaseName_ || fieldName == "alphaRho." + sourcePhaseName_)
    {
        eqn -= mDot;
    }
    else
    {
        for (const targetPhase& tp : targetPhases_)
        {
            if (fieldName == "rho." + tp.name || fieldName == "alphaRho." + tp.name)
            {
                eqn += mDot * tp.fraction;
                return;
            }
            
            size_t dotIndex = fieldName.find('.');
            if (dotIndex != string::npos)
            {
                word phaseName = fieldName.substr(dotIndex + 1);
                if (phaseName == tp.name)
                {
                    word specieName = fieldName.substr(0, dotIndex);
                    if (tp.targetFractions.found(specieName))
                    {
                        eqn += mDot * tp.fraction * tp.targetFractions[specieName];
                        return;
                    }
                }
            }
        }
    }
}


void pressurePhaseChange::addSup
(
    const volVectorField& field,
    fvMatrix<vector>& eqn
) const
{
    // No momentum source implemented yet
}


void pressurePhaseChange::addSup
(
    const volScalarField& rho,
    const volScalarField& field,
    fvMatrix<scalar>& eqn
) const
{
    if (mesh().foundObject<fvTimeIntegrator>("timeIntegrator"))
    {
        return;
    }

    convertMoleFractions();

    const word& fieldName = field.name();
    const volScalarField& alphaSource = 
        mesh().lookupObject<volScalarField>("alpha." + sourcePhaseName_);
    const volScalarField& rhoSource = 
        mesh().lookupObject<volScalarField>("rho." + sourcePhaseName_);

    tmp<volScalarField> tgamma = gamma();
    
    volScalarField mDot = tgamma() * alphaSource * rhoSource / tau_;
    dimensionedScalar residualRho("residualRho", dimDensity, 1e-10);

    if (fieldName == "rhoE" || fieldName == "e" || fieldName == "he")
    {
        tmp<volScalarField> tEOff = energyOffset();
        const volScalarField& eOff = tEOff();
        eqn -= mDot * eOff;
    }
    else if (fieldName == "alpha." + sourcePhaseName_)
    {
        eqn -= mDot / max(rhoSource, residualRho);
    }
    else
    {
        for (const targetPhase& tp : targetPhases_)
        {
            if (fieldName == "alpha." + tp.name)
            {
                eqn += mDot * tp.fraction / max(rhoSource, residualRho);
                break;
            }
        }
    }
}


bool pressurePhaseChange::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        coeffs().lookup("sourcePhase") >> sourcePhaseName_;
        
        targetPhases_.clear();
        if (coeffs().found("targets"))
        {
            const dictionary& targetsDict = coeffs().subDict("targets");
            targetPhases_.setSize(targetsDict.size());
            label i = 0;
            forAllConstIter(dictionary, targetsDict, iter)
            {
                if (iter().isDict())
                {
                    const dictionary& targetDict = iter().dict();
                    targetPhases_[i].name = iter().keyword();
                    targetPhases_[i].fraction = targetDict.lookupOrDefault<scalar>("fraction", 1.0);
                    
                    targetPhases_[i].targetFractions.clear();
                    if (targetDict.found("targetFractions"))
                    {
                        const dictionary& fractionsDict = targetDict.subDict("targetFractions");
                        forAllConstIter(dictionary, fractionsDict, fIter)
                        {
                            targetPhases_[i].targetFractions.insert(fIter().keyword(), readScalar(fIter().stream()));
                        }
                    }
                    i++;
                }
            }
            targetPhases_.setSize(i);
        }
        else
        {
            targetPhases_.setSize(1);
            coeffs().lookup("targetPhase") >> targetPhases_[0].name;
            targetPhases_[0].fraction = 1.0;
            
            targetPhases_[0].targetFractions.clear();
            if (coeffs().found("targetFractions"))
            {
                const dictionary& fractionsDict = coeffs().subDict("targetFractions");
                forAllConstIter(dictionary, fractionsDict, iter)
                {
                    targetPhases_[0].targetFractions.insert(iter().keyword(), readScalar(iter().stream()));
                }
            }
        }

        // Verify that target phase fractions sum to 1.0 (within a small tolerance)
        scalar sumFractions = 0.0;
        for (const targetPhase& tp : targetPhases_)
        {
            sumFractions += tp.fraction;
        }
        if (mag(sumFractions - 1.0) > 1e-5)
        {
            WarningInFunction
                << "Sum of target phase fractions = " << sumFractions
                << ", which is not 1.0. This might lead to mass conservation issues."
                << endl;
        }

        // Configure pressure interval
        pLower_ = 0.0;
        pUpper_ = 1e7;
        pRange_ = 1e6;

        if (coeffs().found("pInterval"))
        {
            scalarList pInt(coeffs().lookup("pInterval"));
            if (pInt.size() == 2)
            {
                pLower_ = min(pInt[0], pInt[1]);
                pUpper_ = max(pInt[0], pInt[1]);
                pRange_ = max(0.05 * (pUpper_ - pLower_), 1.0);
            }
        }
        else if (coeffs().found("pMin") && coeffs().found("pMax"))
        {
            coeffs().lookup("pMin") >> pLower_;
            coeffs().lookup("pMax") >> pUpper_;
            pRange_ = max(0.05 * (pUpper_ - pLower_), 1.0);
        }
        else if (coeffs().found("pLower") && coeffs().found("pUpper"))
        {
            coeffs().lookup("pLower") >> pLower_;
            coeffs().lookup("pUpper") >> pUpper_;
            pRange_ = max(0.05 * (pUpper_ - pLower_), 1.0);
        }
        else if (coeffs().found("pThreshold"))
        {
            coeffs().lookup("pThreshold") >> pUpper_;
            pLower_ = coeffs().lookupOrDefault<scalar>("pLower", 0.0);
            pRange_ = max(0.05 * (pUpper_ - pLower_), 1.0);
        }
        else if (coeffs().found("Tthreshold"))
        {
            // Legacy fallback if an old dictionary specifies Tthreshold
            WarningInFunction
                << "Found legacy 'Tthreshold' in dictionary for " << name()
                << ". Phase change is now pressure-driven ('pressurePhaseChange'). "
                << "Defaulting pressure interval to [0, 100 bar]." << endl;
            pUpper_ = 1e7;
            pLower_ = 0.0;
            pRange_ = 1e6;
        }

        if (coeffs().found("pRange"))
        {
            coeffs().lookup("pRange") >> pRange_;
        }
        else if (coeffs().found("Prange"))
        {
            coeffs().lookup("Prange") >> pRange_;
        }

        if (pRange_ <= 0.0)
        {
            pRange_ = max(1e-4 * pUpper_, 1.0);
        }

        tau_ = coeffs().lookupOrDefault<scalar>("tau", 1e-4);
        isMoleFraction_ =
            coeffs().lookupOrDefault<bool>("isMoleFraction", false);
        fractionsConverted_ = false;
        scalar eOffsetVal = coeffs().lookupOrDefault<scalar>("eOffset", 0.0);
        eOffset_ = dimensionedScalar("eOffset", dimEnergy/dimMass, eOffsetVal);

        Info<< "    " << name() << " (pressurePhaseChange):" << nl
            << "        sourcePhase: " << sourcePhaseName_ << nl
            << "        pInterval: [" << pLower_/1e5 << ", " << pUpper_/1e5 << "] bar (pRange = " << pRange_/1e5 << " bar)" << nl
            << "        tau: " << tau_ << " s" << endl;

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// ************************************************************************* //
