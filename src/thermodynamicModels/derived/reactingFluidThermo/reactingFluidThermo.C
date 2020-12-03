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

#include "reactingFluidThermo.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::reactingFluidThermo<ThermoType>::reactingFluidThermo
(
    const word& name,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
:
    speciesTable(),
    autoPtr<chemistryReader<ThermoType>>
    (
        chemistryReader<ThermoType>::New
        (
            dict,
            *this
        )
    ),
    multicomponentFluidThermo<ThermoType>
    (
        autoPtr<chemistryReader<ThermoType>>::operator()().speciesThermo(),
        name,
        p,
        rho,
        e,
        T,
        dict,
        master,
        masterName
    ),
    PtrList<Reaction<ThermoType>>
    (
        autoPtr<chemistryReader<ThermoType>>::operator()().reactions()
    ),
    speciesComposition_
    (
        autoPtr<chemistryReader<ThermoType>>::operator()().specieComposition()
    ),
    chemistryPtr_(BasicChemistryModel<multicomponentFluidThermoModel>::New(*this))
{
    autoPtr<chemistryReader<ThermoType>>::clear();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::reactingFluidThermo<ThermoType>::~reactingFluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::reactingFluidThermo<ThermoType>::postUpdate()
{
    const fvMesh& mesh(this->rho_.mesh());
    chemistryPtr_->solve(mesh.time().deltaTValue());
    const volScalarField& rho = this->rho();
    const dimensionedScalar& dT(rho.time().deltaT());

    volScalarField YTot(volScalarField::New("YTot", rho.mesh(), 0.0));

    forAll(this->Ys_, i)
    {
        if (this->active(i) && i != this->inertIndex_)
        {
            fvScalarMatrix YEqn
            (
                fvm::ddt(rho, this->Ys_[i])
              - fvc::ddt(rho, this->Ys_[i])
              - fvm::laplacian(this->mu(), this->Ys_[i])
             ==
                chemistryPtr_->RR(i)
            );
            YEqn.solve("Yi");

            this->Ys_[i].max(0.0);
            this->Ys_[i].correctBoundaryConditions();
            YTot += this->Ys_[i];
        }
    }
    this->Ys_[this->inertIndex_] = max(1.0 - YTot, 0.0);
    this->Ys_[this->inertIndex_].correctBoundaryConditions();
    this->e_ += dT*chemistryPtr_->Qdot()/rho;
}

// ************************************************************************* //
