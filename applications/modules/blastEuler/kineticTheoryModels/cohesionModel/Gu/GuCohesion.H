/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
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
    Foam::kineticTheoryModels::cohesionModels::Gu

Description

SourceFiles
    GuCohesion.C

\*---------------------------------------------------------------------------*/

#ifndef GuCohesion_H
#define GuCohesion_H

#include "cohesionModel.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace cohesionModels
{

/*---------------------------------------------------------------------------*\
                           Class Gu Declaration
\*---------------------------------------------------------------------------*/

class Gu
:
    public cohesionModel
{
// Private member data

    //- Maximum Cohesive force
    dimensionedScalar Fmax_;

    //- Minimum volume fraction for cohesive forces
    dimensionedScalar alphaMinCohesion_;

    //- Limit for alpha
    dimensionedScalar alphaMax_;

    //- Coefficients
    dimensionedScalar aCoh1_;
    dimensionedScalar aCoh2_;
    dimensionedScalar aCoh3_;
    dimensionedScalar ae_;
    dimensionedScalar aw_;

    //- Yield stress
    volScalarField tauYield_;

    //- Lambert-W function
    volScalarField W_;

    //- polynomial
    static const Polynomial<5> fPhi_;




public:

    //- Runtime type information
    TypeName("Gu");


    // Constructors

        //- Construct from components
        Gu
        (
            const dictionary& dict,
            const kineticTheoryModel& kt
        );


    //- Destructor
    virtual ~Gu();


     // Member Functions

        //- Update fields
        virtual void update();

        //- Viscosity
        virtual tmp<volScalarField> nu() const;

        //- Granular pressure
        virtual tmp<volScalarField> Ps() const;

        //- Derivative of the granular pressure w.r.t. volume  fraction
        virtual tmp<volScalarField> dPsdAlpha() const;

        //- Derivative of the granular pressure w.r.t. granular temperature
        virtual tmp<volScalarField> dPsdTheta() const;

        //- Dissipation source
        virtual tmp<volScalarField> dissipationSource
        (
            const dimensionedScalar&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cohesionModels
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
