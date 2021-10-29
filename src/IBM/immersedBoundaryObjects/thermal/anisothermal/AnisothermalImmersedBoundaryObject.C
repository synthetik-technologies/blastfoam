/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
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

#include "AnisothermalImmersedBoundaryObject.H"
#include "fvm.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::
AnisothermalImmersedBoundaryObject
(
    const polyMesh& mesh,
    const dictionary& dict,
    const dictionary& stateDict
)
:
    ImmersedType(mesh, dict, stateDict),
    thermo_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::
~AnisothermalImmersedBoundaryObject()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ImmersedType>
void Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::initialize()
{
    ImmersedType::initialize();
    const fvMesh& mesh = this->immersedMesh();
    if (mesh.foundObject<solidBlastThermo>(basicThermo::dictName))
    {
        thermo_.set
        (
            &mesh.lookupObjectRef<solidBlastThermo>(basicThermo::dictName)
        );
    }
    else
    {
        thermo_.set
        (
            solidBlastThermo::New
            (
                mesh,
                IOdictionary
                (
                    IOobject
                    (
                        basicThermo::dictName,
                        mesh.time().constant(),
                        this->immersedMesh(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                )
            ).ptr()
        );
    }
}

template<class ImmersedType>
void Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::solve()
{
    if (!ImmersedType::solveTemperature())
    {
        const volScalarField& rho(thermo_->rho());
        Foam::solve
        (
            fvm::ddt(rho, thermo_->Cp(), thermo_->T())
          ==
            fvm::laplacian(thermo_->kappa(), thermo_->T())
        );
    }

    ImmersedType::solve();
}


template<class ImmersedType>
void Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::status() const
{
    ImmersedType::status();
    const volScalarField& T(thermo_->T());

    Info<< "    Temperature (min, max, avg): "
        << min(T).value() << ", "
        << max(T).value() << ", "
        << T.weightedAverage(this->immersedMesh().V()).value()
        << endl;
}


template<class ImmersedType>
void Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::movePoints()
{
    ImmersedType::movePoints();
}


template<class ImmersedType>
bool Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::read
(
    const dictionary& dict
)
{
    return ImmersedType::read(dict);
}


template<class ImmersedType>
void Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::write
(
    Ostream& os
) const
{
    ImmersedType::write(os);
}


template<class ImmersedType>
void Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::write
(
    dictionary& dict
) const
{
    ImmersedType::write(dict);
}

// ************************************************************************* //
