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


tmp<volScalarField> pressurePhaseChange::energyOffset() const
{
    tmp<volScalarField> teOffset
    (
        volScalarField::New
        (
            "energyOffset",
            mesh(),
            eOffset_
        )
    );
    volScalarField& eOffset = teOffset.ref();

    if (mesh().foundObject<multiphaseFluidBlastThermo>("thermophysicalProperties"))
    {
        const multiphaseFluidBlastThermo& thermo = 
            mesh().lookupObject<multiphaseFluidBlastThermo>("thermophysicalProperties");

        const wordList& phaseNames = thermo.phaseNames();
        
        label sourcePhasei = -1;
        forAll(phaseNames, i)
        {
            if (phaseNames[i] == sourcePhaseName_)
            {
                sourcePhasei = i;
                break;
            }
        }

        if (sourcePhasei != -1)
        {
            const volScalarField& T = mesh().lookupObject<volScalarField>("T");
            scalarField& eOffsetCells = eOffset.primitiveFieldRef();
            
            forAll(eOffsetCells, celli)
            {
                scalar TCell = T[celli];
                scalar eSource = thermo.thermo(sourcePhasei).cellhe(TCell, celli);
                scalar eTargets = 0.0;
                for (const targetPhase& tp : targetPhases_)
                {
                    label targetPhasei = -1;
                    forAll(phaseNames, i)
                    {
                        if (phaseNames[i] == tp.name)
                        {
                            targetPhasei = i;
                            break;
                        }
                    }
                    if (targetPhasei != -1)
                    {
                        eTargets += tp.fraction * thermo.thermo(targetPhasei).cellhe(TCell, celli);
                    }
                }
                eOffsetCells[celli] = eTargets - eSource;
            }
            
            forAll(eOffset.boundaryField(), patchi)
            {
                scalarField& peOffset = eOffset.boundaryFieldRef()[patchi];
                const scalarField& pT = T.boundaryField()[patchi];
                forAll(peOffset, facei)
                {
                    scalar TFace = pT[facei];
                    scalar eSource = thermo.thermo(sourcePhasei).patchFacehe(TFace, patchi, facei);
                    scalar eTargets = 0.0;
                    for (const targetPhase& tp : targetPhases_)
                    {
                        label targetPhasei = -1;
                        forAll(phaseNames, i)
                        {
                            if (phaseNames[i] == tp.name)
                            {
                                targetPhasei = i;
                                break;
                            }
                        }
                        if (targetPhasei != -1)
                        {
                            eTargets += tp.fraction * thermo.thermo(targetPhasei).patchFacehe(TFace, patchi, facei);
                        }
                    }
                    peOffset[facei] = eTargets - eSource;
                }
            }
        }
    }

    return teOffset;
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

        // Energy source for rhoE
        tmp<volScalarField> tEOff = energyOffset();
        const volScalarField& eOff = tEOff();
        volScalarField::Internal SrhoE
        (
            IOobject("SrhoE", mesh().time().name(), mesh()),
            mesh(),
            dimensionedScalar("SrhoE", dimEnergy/dimVolume/dimTime, 0.0)
        );
        SrhoE.primitiveFieldRef() = -mDot.primitiveField() * eOff.primitiveField();
        integrator.addSource("rhoE", SrhoE);

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
        coeffs().lookupOrDefault<bool>("isMoleFraction", isMoleFraction_);
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
