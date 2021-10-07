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
    powerLawPlastic

Description
    Power law plastic mechanical law:

        \sigma_y = k \epsilon^n

SourceFiles
    powerLawPlastic.C


\*---------------------------------------------------------------------------*/

#ifndef powerLawPlastic_H
#define powerLawPlastic_H

#include "plasticModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class powerLawPlastic Declaration
\*---------------------------------------------------------------------------*/

template<class PlasticType>
class powerLawPlastic
:
    public PlasticType
{
    // Private data

        //- Strength coefficient
        dimensionedScalar k_;

        //- Strain hardening exponentn
        dimensionedScalar n_;

        //- Strain at yield
        dimensionedScalar epsilonY_;

    // Private Member Functions

        //- Return the current Kirchhoff yield stress
        virtual scalar curYieldStress
        (
            const scalar curEpsilonPEq,     // Current equivalent plastic strain
            const scalar curEpsilonPEqDot,  // Current equivalent plastic strain rate
            const scalar J = 1.0            // Current Jacobian
        ) const;

        virtual bool nonLinearPlasticity() const
        {
            return true;
        }

public:


    // Constructors

        //- Construct from dictionary
        powerLawPlastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~powerLawPlastic();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "powerLawPlastic.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
