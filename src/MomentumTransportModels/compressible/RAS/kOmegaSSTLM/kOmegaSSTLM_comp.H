/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2025-06-09 Jeff Heylmun     : Derived from compressible::kOmegaSST
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

Class
    Foam::RASModels::kOmegaSSTLM_comp

Description
    Langtry-Menter 4-equation transitional SST model
    based on the k-omega-SST RAS model.

    References:
    \verbatim
        Langtry, R. B., & Menter, F. R. (2009).
        Correlation-based transition modeling for unstructured parallelized
        computational fluid dynamics codes.
        AIAA journal, 47(12), 2894-2906.

        Menter, F. R., Langtry, R., & Volker, S. (2006).
        Transition modelling for general purpose CFD codes.
        Flow, turbulence and combustion, 77(1-4), 277-303.

        Langtry, R. B. (2006).
        A correlation-based transition model using local variables for
        unstructured parallelized CFD codes.
        Phd. Thesis, Universität Stuttgart.
    \endverbatim

    The model coefficients are
    \verbatim
        compressible::kOmegaSSTLMCoeffs
        {
            // Default SST coefficients
            alphaK1     0.85;
            alphaK2     1;
            alphaOmega1 0.5;
            alphaOmega2 0.856;
            beta1       0.075;
            beta2       0.0828;
            betaStar    0.09;
            gamma1      5/9;
            gamma2      0.44;
            a1          0.31;
            b1          1;
            c1          10;
            F3          no;

            // Default LM coefficients
            ca1         2;
            ca2         0.06;
            ce1         1;
            ce2         50;
            cThetat     0.03;
            sigmaThetat 2;

            lambdaErr   1e-6;
            maxLambdaIter 10;
        }
    \endverbatim

SourceFiles
    kOmegaSSTLM_comp.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTLM_comp_H
#define kOmegaSSTLM_comp_H

#include "kOmegaSST_comp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                         Class kOmegaSSTLM_comp Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kOmegaSSTLM_comp
:
    public kOmegaSST_comp<BasicMomentumTransportModel>
{
protected:

    // Protected data

        // Model constants

            dimensionedScalar ca1_;
            dimensionedScalar ca2_;

            dimensionedScalar ce1_;
            dimensionedScalar ce2_;

            dimensionedScalar cThetat_;
            dimensionedScalar sigmaThetat_;

            //- Convergence criterion for the lambda/thetat loop
            scalar lambdaErr_;

            //- Maximum number of iterations to converge the lambda/thetat loop
            label maxLambdaIter_;

            //- Stabilisation for division by the magnitude of the velocity
            const dimensionedScalar stabU_;

            //- Stabilisation for division by the boundary layer thickness
            const dimensionedScalar stabDelta_;


        // Fields

            //- Transition onset momentum-thickness Reynolds number
            volScalarField ReThetat_;

            //- Intermittency
            volScalarField gammaInt_;

            //- Effective intermittency
            volScalarField::Internal gammaIntEff_;

            //- Transition onset momentum-thickness Reynolds number
            // (based on freestream conditions)
            mutable volScalarField::Internal lambda_;


    // Protected Member Functions

        //- Modified form of the k-omega SST F1 function
        virtual tmp<volScalarField> F1(const volScalarField& CDkOmega) const;

        //- Modified form of the k-omega SST k production rate
        virtual tmp<volScalarField::Internal> Pk
        (
            const volScalarField::Internal& G
        ) const;

        //- Modified form of the k-omega SST epsilon/k
        virtual tmp<volScalarField::Internal> epsilonByk
        (
            const volScalarField::Internal& F1,
            const volScalarField::Internal& F2
        ) const;

        //- Freestream blending-function
        tmp<volScalarField::Internal> Fthetat
        (
            const volScalarField::Internal& Us,
            const volScalarField::Internal& Omega,
            const volScalarField::Internal& nu
        ) const;

        //- Empirical correlation for critical Reynolds number where the
        //  intermittency first starts to increase in the boundary layer
        tmp<volScalarField::Internal> ReThetac() const;

        //- Empirical correlation that controls the length of the
        //  transition region
        tmp<volScalarField::Internal> Flength
        (
            const volScalarField::Internal& nu
        ) const;

        //- Transition onset location control function
        tmp<volScalarField::Internal> Fonset
        (
            const volScalarField::Internal& Rev,
            const volScalarField::Internal& ReThetac,
            const volScalarField::Internal& RT
        ) const;

        //- Return the transition onset momentum-thickness Reynolds number
        // (based on freestream conditions)
        tmp<volScalarField::Internal> ReThetat0
        (
            const volScalarField::Internal& Us,
            const volScalarField::Internal& dUsds,
            const volScalarField::Internal& nu
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        void correctReThetatGammaInt();


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    CompressibleTypeName("kOmegaSSTLM");


    // Constructors

        //- Construct from components
        kOmegaSSTLM_comp
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = baseName()
        );

        //- Disallow default bitwise copy construction
        kOmegaSSTLM_comp(const kOmegaSSTLM_comp&) = delete;


    //- Destructor
    virtual ~kOmegaSSTLM_comp()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Access function transition onset momentum-thickness Reynolds number
        const volScalarField& ReThetat() const
        {
            return ReThetat_;
        }

        //- Access function to intermittency
        const volScalarField& gammaInt() const
        {
            return gammaInt_;
        }

        //- Return the effective diffusivity for transition onset
        //  momentum-thickness Reynolds number
        tmp<volScalarField> DReThetatEff() const
        {
            return volScalarField::New
            (
                "DReThetatEff",
                sigmaThetat_*(this->nut_ + this->nu())
            );
        }

        //- Return the effective diffusivity for intermittency
        tmp<volScalarField> DgammaIntEff() const
        {
            return volScalarField::New
            (
                "DgammaIntEff",
                this->nut_ + this->nu()
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const kOmegaSSTLM_comp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kOmegaSSTLM_comp.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
