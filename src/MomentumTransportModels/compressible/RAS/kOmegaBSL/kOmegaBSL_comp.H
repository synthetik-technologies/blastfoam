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
    Foam::RASModels::kOmegaBSL_comp

Description
    Base line k-omega model

    References:
    \verbatim
        //TODO
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

            lambdaErr   1e-6;
            maxLambdaIter 10;
        }
    \endverbatim

SourceFiles
    kOmegaBSL_comp.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaBSL_comp_H
#define kOmegaBSL_comp_H

#include "kOmegaSSTBase_comp.H"
#include "eddyViscosity.H"
#include "RASModel.H"
#include "compressibilityCorrection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                         Class kOmegaBSL_comp Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kOmegaBSL_comp
:
    public ::Foam::compressible::kOmegaSST
    <
        eddyViscosity<RASModel<BasicMomentumTransportModel>>,
        BasicMomentumTransportModel
    >
{
protected:

        virtual void correctNut
        (
            const volScalarField& S2,
            const volScalarField& F2
        );
        virtual void correctNut();

        // //- Return k production rate
        // virtual tmp<volScalarField::Internal> Pk
        // (
        //     const volScalarField::Internal& G
        // ) const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    CompressibleTypeName("kOmegaBSL");


    // Constructors

        //- Construct from components
        kOmegaBSL_comp
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = baseName(),
            const HashTable<scalar>& defaults = HashTable<scalar>()
        );

        //- Disallow default bitwise copy construction
        kOmegaBSL_comp(const kOmegaBSL_comp&) = delete;


    //- Destructor
    virtual ~kOmegaBSL_comp()
    {}


    // Member Functions

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Disallow default bitwise assignment
        void operator=(const kOmegaBSL_comp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kOmegaBSL_comp.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
