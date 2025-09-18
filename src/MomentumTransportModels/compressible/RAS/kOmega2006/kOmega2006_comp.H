/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

Class
    Foam::RASModels::kOmega2006_comp

Description
    Standard (2006) high Reynolds-number k-omega turbulence model for
    incompressible and compressible flows.

    References:
    \verbatim
        Wilcox, D. C. (2006).
        Turbulence modeling for CFD, 3rd edition
        La Canada, CA: DCW industries, Inc.

        Wilcox, D. C. (2008).
        Formulation of the kw turbulence model revisited.
        AIAA journal, 46(11), 2823-2838.
    \endverbatim

    The default model coefficients are
    \verbatim
        compressible::kOmega2006Coeffs
        {
            Cmu         0.09;
            beta0       0.0708;
            gamma       0.52;
            Clim        0.875;
            alphak      0.6;
            alphaOmega  0.5;
        }
    \endverbatim

SourceFiles
    kOmega2006_comp.C

\*---------------------------------------------------------------------------*/

#ifndef kOmega2006_comp_H
#define kOmega2006_comp_H

#include "kOmega2006.H"
#include "compressibilityCorrection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kOmega2006_comp Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kOmega2006_comp
:
    public kOmega2006<BasicMomentumTransportModel>,
    public ::Foam::compressible::correction
{
    // Private Member Functions

        void correctNut(const volTensorField& gradU);

        tmp<volScalarField::Internal> beta(const volTensorField& gradU) const;

        tmp<volScalarField::Internal> CDkOmega() const;

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("compressible::kOmega2006");


    // Constructors

        //- Construct from components
        kOmega2006_comp
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = "kOmega2006"
        );


    //- Destructor
    virtual ~kOmega2006_comp()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmega2006_comp.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
