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
    Foam::equationOfStates::Tait

Description

    References:
    \verbatim
        Zheng, H.W., Shu, C., Chew, Y.T., Qin, N.  (2011).
        A solution adaptive simulation of compressible multi-fluid flows
        with general equation of state.
        International Journal for Numerical Methods in Fluids, 67, 616-637.
    \endverbatim

SourceFiles
    Tait.C

\*---------------------------------------------------------------------------*/

#ifndef Tait_H
#define Tait_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationOfState.H"

namespace Foam
{
namespace equationOfStates
{

/*---------------------------------------------------------------------------*\
                           Class Tait Declaration
\*---------------------------------------------------------------------------*/

class Tait
:
    public equationOfState
{
// Private data

        //- Model coefficients
        dimensionedScalar a_;
        dimensionedScalar b_;

        //- Specific heat ratio
        dimensionedScalar gamma_;



public:

    //- Runtime type information
    ClassName("Tait");

    // Constructor
    Tait(const volScalarField& e, const dictionary& dict);


    //- Destructor
    virtual ~Tait();


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
