/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    Foam::saturationPressureModels::AntoineExtended

Description
    Constant saturation pressure

SourceFiles
    AntoineExtendedSaturationPressureModel.C

\*---------------------------------------------------------------------------*/

#ifndef AntoineExtendedSaturationPressureModel_H
#define AntoineExtendedSaturationPressureModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "saturationPressureModel.H"

namespace Foam
{
namespace saturationPressureModels
{

/*---------------------------------------------------------------------------*\
                           Class AntoineExtended Declaration
\*---------------------------------------------------------------------------*/

class AntoineExtended
:
    public saturationPressureModel
{
    // Private data

        //- Model coefficients
        dimensionedScalar A_;
        dimensionedScalar B_;
        dimensionedScalar C_;
        dimensionedScalar D_;
        dimensionedScalar E_;
        dimensionedScalar F_;


public:

    //- Runtime type information
    TypeName("AntoineExtended");

    // Constructor
    AntoineExtended(const dictionary& dict);


    //- Destructor
    virtual ~AntoineExtended();


    // Member Functions

        //- Return satuation pressure given a temperature
        virtual scalar pSat(const scalar T) const;

        //- Return derivative of satuation pressure w.r.t temperature
        virtual scalar derivative(const scalar T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationPressureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
