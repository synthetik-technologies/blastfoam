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
    Foam::saturationTemperatureModels::constantTemperature

Description
    Constant saturation temperature model.

SourceFiles
    constantSaturationTemperatureModel.C

\*---------------------------------------------------------------------------*/

#ifndef constantSaturationTemperatureModel_H
#define constantSaturationTemperatureModel_H

#include "saturationTemperatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationTemperatureModels
{

/*---------------------------------------------------------------------------*\
                        Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public saturationTemperatureModel
{
protected:

    // Private Data

        //- Constant saturation temperature
        const dimensionedScalar Tsat_;

public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from a dictionary
        constant(const dictionary& dict);

        //- Construct from a value
        constant(const dimensionedScalar& Tsat);


    //- Destructor
    virtual ~constant();


    // Member Functions

        using saturationTemperatureModel::Tsat;

        //- Saturation temperature
        virtual scalar Tsat(const scalar p) const
        {
            return Tsat_.value();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationTemperatureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
