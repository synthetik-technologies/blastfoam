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
    Foam::idealGas

Description

SourceFiles
    idealGas.C

\*---------------------------------------------------------------------------*/

#ifndef idealGas_H
#define idealGas_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationOfState.H"

namespace Foam
{
namespace equationOfStates
{

/*---------------------------------------------------------------------------*\
                           Class idealGas Declaration
\*---------------------------------------------------------------------------*/

class idealGas
:
    public equationOfState
{
// Private data

        //- Specific heat ratio
        dimensionedScalar gamma_;



public:

    //- Runtime type information
    ClassName("idealGas");

    // Constructor
    idealGas(const volScalarField& e, const dictionary& dict);


    //- Destructor
    virtual ~idealGas();


    // Member Functions

        //- Gamma function
        virtual tmp<volScalarField> Gamma() const;

        //- Pi function
        virtual tmp<volScalarField> Pi() const;

        //- Delta coefficient
        virtual tmp<volScalarField> delta() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace equationOfStates
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
