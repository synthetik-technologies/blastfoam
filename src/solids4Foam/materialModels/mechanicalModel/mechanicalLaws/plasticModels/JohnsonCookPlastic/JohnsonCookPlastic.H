/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    JohnsonCookPlastic

Description
    Johnson-Cook plastic mechanical law:

        tau_ = (A + B*\bar{\epsilon}^{p*})*(1 + C*ln(\dot{\epsilon}^{*}))(1-T^{*n})

SourceFiles
    JohnsonCookPlastic.C


\*---------------------------------------------------------------------------*/

#ifndef JohnsonCookPlastic_H
#define JohnsonCookPlastic_H

#include "plasticModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class JohnsonCookPlastic Declaration
\*---------------------------------------------------------------------------*/

template<class PlasticType>
class JohnsonCookPlastic
:
    public PlasticType
{
    // Private data

        //- Const reference to temperature
        const volScalarField& T_;

        //- Saved temporary interpolated temperature
        tmp<surfaceScalarField> Tf_;

        //- Saved effective elastic strain rate
        tmp<volScalarField> epsilonElasticEffDot_;
        tmp<surfaceScalarField> epsilonElasticEffDotf_;

        //- Yield stress at reference conditions
        dimensionedScalar A_;

        //- Strain hardening constant
        dimensionedScalar B_;

        //- Strengthing coefficient of strain rate
        dimensionedScalar C_;

        //- Strain hardening exponentn
        dimensionedScalar n_;

        //- Thermal softening coefficient
        dimensionedScalar m_;

        //- Reference temperature
        dimensionedScalar Tref_;

        //- Melting temperature
        dimensionedScalar Tm_;

        //- Reference plastic strain
        dimensionedScalar epsilonP0_;

        //- Referenece plastic strain rate
        dimensionedScalar epsilonPDot0_;

        //- Is the total strain rate used
        Switch totalStrainRate_;

        //- Are temperature effects included
        Switch temperatureEffects_;

        mutable scalar Ti_;
        mutable scalar eDoti_;

    // Private Member Functions

        //- Return the current Kirchhoff yield stress
        virtual scalar curYieldStress
        (
            const scalar curEpsilonPEq,     // Current equivalent plastic strain
            const scalar curEpsilonPEqDot,  // Current equivalent plastic strain rate
            const scalar J = 1.0            // Current Jacobian
        ) const;

        //- Set local values for fields not passed in functions
        virtual void setCellValues(const label celli) const;
        virtual void setVolPatchFaceValues
        (
            const label patchi,
            const label facei
        ) const;
        virtual void setFaceValues(const label facei) const;
        virtual void setSurfacePatchFaceValues
        (
            const label patchi,
            const label facei
        ) const;

        virtual bool nonLinearPlasticity() const
        {
            return true;
        }

public:


    // Constructors

        //- Construct from dictionary
        JohnsonCookPlastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~JohnsonCookPlastic();


    // Member Functions

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "JohnsonCookPlastic.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
