/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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

#include "reactingFluidBlastThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::reactingFluidBlastThermo<Thermo>::reactingFluidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    multicomponentFluidBlastThermo<Thermo>
    (
        mesh,
        dict,
        phaseName,
        masterName
    ),
    combustion_(blastCombustionModel::NewFluid(*this)),
    odeCombustion_
    (
        dict.lookupOrDefaultBackwardsCompatible
        (
            {"odeCombustion", "odeChemistry"},
            false
        )
    ),
    Qdot_
    (
        IOobject
        (
            IOobject::groupName("Qdot", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    )
{
    this->correctMassFractions();
}


template<class Thermo>
Foam::reactingFluidBlastThermo<Thermo>::reactingFluidBlastThermo
(
    const HashPtrTable<typename Thermo::thermoType, word, string::hash>& thermoData,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    multicomponentFluidBlastThermo<Thermo>
    (
        thermoData,
        mesh,
        dict,
        phaseName,
        masterName
    ),
    combustion_(blastCombustionModel::NewFluid(*this)),
    odeCombustion_
    (
        dict.lookupOrDefaultBackwardsCompatible
        (
            {"odeCombustion", "odeChemistry"},
            false
        )
    ),
    Qdot_
    (
        IOobject
        (
            IOobject::groupName("Qdot", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    )
{
    this->correctMassFractions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::reactingFluidBlastThermo<Thermo>::~reactingFluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::reactingFluidBlastThermo<Thermo>::update()
{
    if (this->step() == 1 || odeCombustion_)
    {
        combustion_->correct();
    }

    Qdot_ = combustion_->Qdot();

    forAll(this->species_, speciei)
    {
        if (this->solveSpecie(speciei))
        {
            this->addDelta(this->species_[speciei], combustion_->R(speciei));
        }
    }
    multicomponentFluidBlastThermo<Thermo>::update();
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::reactingFluidBlastThermo<Thermo>::ESource() const
{
    return Qdot_;
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::reactingFluidBlastThermo<Thermo>::ESource
(
    const volScalarField& alpha
) const
{
    return alpha*ESource();
}

// ************************************************************************* //
