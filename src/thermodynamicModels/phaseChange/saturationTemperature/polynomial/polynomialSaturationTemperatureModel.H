/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
28-08-2025  Synthetik Applied Technologies: Added single input functions
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
    Foam::saturationTemperatureModels::polynomial

Description
    Polynomial equation for the saturation vapour temperature in terms of
    the vapour pressure (in Pa).

    \f[
        T_sat = \sum_i C_i p^i
    \f]

    where \f$p\f$ is the pressure in Pa and \f$C\f$ are the coefficients.

SourceFiles
    polynomialSaturationTemperatureModel.C

\*---------------------------------------------------------------------------*/

#ifndef polynomialSaturationTemperatureModel_H
#define polynomialSaturationTemperatureModel_H

#include "saturationTemperatureModel.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationTemperatureModels
{

/*---------------------------------------------------------------------------*\
                         Class polynomial Declaration
\*---------------------------------------------------------------------------*/

class polynomial
:
    public saturationTemperatureModel
{
    // Private Data

        //- Polynomial coefficients
        Polynomial<8> C_;


public:

    //- Runtime type information
    TypeName("polynomial");


    // Constructors

        //- Construct from a dictionary
        polynomial(const dictionary& dict);


    //- Destructor
    virtual ~polynomial();


    // Member Functions

        using saturationTemperatureModel::Tsat;

        //- Saturation temperature for volScalarField::Internal
        virtual scalar Tsat(const scalar p) const
        {
            return C_.value(p);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationTemperatureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
