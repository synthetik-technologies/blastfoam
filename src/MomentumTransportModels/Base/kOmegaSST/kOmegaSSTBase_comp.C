/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2025-06-09 Jeff Heylmun     : Derived from compressibilityCorrection
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

#include "kOmegaSSTBase_comp.H"
#include "fluidThermo.H"
#include "fluidBlastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return
        correction::Fcorr(this->k_)
       *::Foam::kOmegaSST
        <
            MomentumTransportModel,
            BasicMomentumTransportModel
        >::Pk(G);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::kOmegaSST
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
:
    ::Foam::kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),
    correction(this->coeffDict_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
bool kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::read()
{
    if (MomentumTransportModel::read())
    {
        return correction::read(this->coeffDict_);
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
