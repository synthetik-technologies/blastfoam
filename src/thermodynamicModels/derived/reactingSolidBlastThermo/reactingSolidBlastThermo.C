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

#include "reactingSolidBlastThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::reactingSolidBlastThermo<Thermo>::reactingSolidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    multicomponentSolidBlastThermo<Thermo>
    (
        mesh,
        dict,
        phaseName,
        masterName
    ),
    combustion_(blastCombustionModel::NewSolid(*this)),
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
Foam::reactingSolidBlastThermo<Thermo>::reactingSolidBlastThermo
(
    const HashPtrTable<typename Thermo::thermoType, word, string::hash>& thermoData,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    multicomponentSolidBlastThermo<Thermo>
    (
        thermoData,
        mesh,
        dict,
        phaseName,
        masterName
    ),
    combustion_(blastCombustionModel::NewSolid(*this)),
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
Foam::reactingSolidBlastThermo<Thermo>::~reactingSolidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Thermo>
void Foam::reactingSolidBlastThermo<Thermo>::solve()
{
    // Chemistry from the state at the start of the time step (first stage,
    // step 0), applied in all stages
    if (this->step() == 0 || odeCombustion_)
    {
        combustion_->correct();
    }

    Qdot_ = combustion_->Qdot();

    // The reaction rates are per unit volume of the phase; the species
    // equations are for alpha*rho*Y (the phase energy source is alpha*Qdot),
    // so the rates are weighted by the volume fraction (one for a single
    // phase)
    const volScalarField alpha
    (
        this->alphaRho()/max(this->rho(), dimensionedScalar(dimDensity, small))
    );

    forAll(this->species_, speciei)
    {
        if (this->solveSpecie(speciei))
        {
            this->addDelta
            (
                this->species_[speciei],
                tmp<volScalarField::Internal>(alpha()*combustion_->R(speciei))
            );
        }
    }
    multicomponentSolidBlastThermo<Thermo>::solve();
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::reactingSolidBlastThermo<Thermo>::ESource() const
{
    return Qdot_;
}


// ************************************************************************* //
