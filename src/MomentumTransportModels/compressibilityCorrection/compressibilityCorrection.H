/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::compressible::correction

Description
    Implementation of the Sarkar/Zemen compressibility correction

    \verbatim
        Sarkar, S., Erlebacher, G., Hussaini, M. Y., and Kreiss, H. O. (2001).
        The Analysis and Modeling of Dilatational Terms in Compressible Turbulence.
        NASA CR 181959
    \endverbatim

    \verbatim
         Zeman, O. (1990).
        Dilatational Dissipation: The Concept and Application in Modeling
        Compressible Mixing Layers.
        Physics of Fluids A, Vol. 2, No. 2,
    \endverbatim

SourceFiles
    compressibilityCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef compressibilityCorrection_H
#define compressibilityCorrection_H

#include "dimensionedScalar.H"
#include "volFields.H"
#include "viscosity.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
                          Class compressible::correction Declaration
\*---------------------------------------------------------------------------*/

class correction
{
public:

    enum MODEL
    {
        K_EPSILON,
        K_OMEGA,
        OTHER
    };

    enum F_CC_Method
    {
        NONE,
        WILCOX,
        ZEMAN,

        SIZE_F
    };
    static const NamedEnum
    <
        F_CC_Method,
        F_CC_Method::SIZE_F
    > F_CC_MethodNames_;

protected:

    // Reference to turbulence model
    const ::Foam::momentumTransportModel& turb_;

    //- Compressibility correction
    F_CC_Method FMethod_;
    dimensionedScalar XiStar_;
    dimensionedScalar lambda_;
    dimensionedScalar Mt0_;


    // Pressure dialation
    bool pressureDilation_;
    dimensionedScalar alphaK2_;
    dimensionedScalar alphaK3_;


public:

    // Constructors

        //- Construct from components
        correction
        (
            dictionary& dict,
            const MODEL model,
            const ::Foam::momentumTransportModel& turb,
            const F_CC_Method ccMethod = F_CC_Method::NONE
        );

        //- Disallow default bitwise copy construction
        correction(const correction&) = delete;


    //- Destructor
    virtual ~correction()
    {}

    // Member Functions

        //- Read from dictionary
        bool read(const dictionary& dict);

        //- Return speed of sound
        tmp<volScalarField::Internal> speedOfSound() const;

        //- Return tubulent Mach number
        tmp<volScalarField::Internal> Mt() const;

        //- Return tubulent Mach number squared
        tmp<volScalarField::Internal> MtSqr() const;

        //- Return tubulent Mach number squared divided by k
        //  for use with implicit sources
        tmp<volScalarField::Internal> MtSqrByk() const;

        //  for use with implicit sources
        tmp<volScalarField::Internal> Fcorr() const;
        tmp<volScalarField::Internal> Fcorr
        (
            tmp<volScalarField::Internal>& tMt
        ) const;
        tmp<volScalarField::Internal> Fcorr
        (
            const volScalarField::Internal& Mt
        ) const;

        //- Set compressible beta and betaStar
        template<class BetaType, class BetaStarType>
        void blendBeta
        (
            const BetaType& beta0,
            const BetaStarType& betaStar0,
            tmp<volScalarField::Internal>& tMt,
            tmp<volScalarField::Internal>& beta,
            tmp<volScalarField::Internal>& betaStar
        ) const
        {
            tmp<volScalarField::Internal> tF(Fcorr(tMt));
            beta = beta0 - betaStar0*XiStar_*tF();
            betaStar = betaStar0*(1.0 + XiStar_*tF());
        }

        //- Pressure dialation source
        tmp<volScalarField::Internal> pressureDialationSource
        (
            const volScalarField::Internal& G
        ) const;
        tmp<volScalarField::Internal> pressureDialationSource
        (
            const volScalarField::Internal& G,
            tmp<volScalarField::Internal>& tMt
        ) const;
};


#define CompressibleTypeName(Name)                                             \
    TypeName("compressible::" Name);                                           \
    static word baseName() {return Name;}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
