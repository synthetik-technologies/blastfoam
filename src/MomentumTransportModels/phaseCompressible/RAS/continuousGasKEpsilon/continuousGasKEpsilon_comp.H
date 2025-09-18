/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Class
    Foam::RASModels::continuousGasKEpsilon_comp

Description
    k-epsilon model for the gas-phase in a two-phase system
    supporting phase-inversion.

    In the limit that the gas-phase fraction approaches zero a contribution from
    the other phase is blended into the k and epsilon equations up to the
    phase-fraction of alphaInversion at which point phase-inversion is
    considered to have occurred and the model reverts to the pure single-phase
    form.

    This model is unpublished and is provided as a stable numerical framework
    on which a more physical model may be built.

    The default model coefficients are
    \verbatim
        continuousGasKEpsilon_compCoeffs
        {
            Cmu             0.09;
            C1              1.44;
            C2              1.92;
            C3              0;
            sigmak          1.0;
            sigmaEps        1.3;
            alphaInversion  0.7;
        }
    \endverbatim

SourceFiles
    continuousGasKEpsilon_comp.C

\*---------------------------------------------------------------------------*/

#ifndef continuousGasKEpsilon_comp_H
#define continuousGasKEpsilon_comp_H

#include "kEpsilon_comp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                    Class continuousGasKEpsilon_comp Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class continuousGasKEpsilon_comp
:
    public kEpsilon_comp<BasicMomentumTransportModel>
{
    // Private Data

        mutable const momentumTransportModel *liquidTurbulencePtr_;

        volScalarField nutEff_;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar alphaInversion_;


    // Protected Member Functions

        virtual void correctNut();
        tmp<volScalarField> phaseTransferCoeff() const;
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("compressible::continuousGasKEpsilon");


    // Constructors

        //- Construct from components
        continuousGasKEpsilon_comp
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = "continuousGasKEpsilon"
        );

        //- Disallow default bitwise copy construction
        continuousGasKEpsilon_comp(const continuousGasKEpsilon_comp&) = delete;


    //- Destructor
    virtual ~continuousGasKEpsilon_comp()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the turbulence model for the liquid phase
        const momentumTransportModel& liquidTurbulence() const;

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const;

        //- Return the effective density for the stress
        virtual tmp<volScalarField> rhoEff() const;

        //- Return the stress tensor [m^2/s^2]
        virtual tmp<volSymmTensorField> sigma() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const continuousGasKEpsilon_comp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "continuousGasKEpsilon_comp.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
