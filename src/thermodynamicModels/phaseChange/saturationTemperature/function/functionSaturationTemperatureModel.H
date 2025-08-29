/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2024 OpenFOAM Foundation
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
    Foam::saturationTemperatureModels::function

Description
    Saturation vapour temperature in terms of
    the vapour pressure (in Pa). The saturation temperature in Kelvins is
    specified as a Foam::Function1 type, to enable use of, e.g. constant,
    polynomial, table values.

    Examples:

    \verbatim
        type function;
        function polynomial
        (
            (308.0422       0)
            (0.0015096      1)
            (-1.61589e-8    2)
            (1.114106e-13   3)
            (-4.52216e-19   4)
            (1.05192e-24    5)
            (-1.2953e-30    6)
            (6.5365e-37     7)
        );
    \endverbatim

    \verbatim
        type function;
        function table;
        functionCoeffs
        {
            file                "filename.csv";
            format              csv;
            nHeaderLine         1;
            columns             (0 1);
            separator           ",";
            mergeSeparators     no;
            outOfBounds         clamp;
            interpolationScheme linear;
        };
    \endverbatim

SourceFiles
    function1SaturationTemperatureModel.C

\*---------------------------------------------------------------------------*/

#ifndef functionSaturationTemperatureModel_H
#define functionSaturationTemperatureModel_H

#include "saturationTemperatureModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationTemperatureModels
{

/*---------------------------------------------------------------------------*\
                          Class function Declaration
\*---------------------------------------------------------------------------*/

class function
:
    public saturationTemperatureModel
{
    // Private Data

        //- Saturation temperature as a function of pressure
        autoPtr<Function1<scalar>> function_;


public:

    //- Runtime type information
    TypeName("function");


    // Constructors

        //- Construct from a dictionary and an interface
        function(const dictionary& dict);


    //- Destructor
    virtual ~function();


    // Member Functions

        using saturationTemperatureModel::Tsat;

        //- Saturation temperature
        virtual scalar Tsat(const scalar p) const
        {
            return function_->value(p);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationTemperatureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
