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
    Foam::equationOfStates::JWL

Description
    Equation of state for an energetic material

    References:
    \verbatim
        Zheng, H.W., Shu, C., Chew, Y.T., Qin, N.  (2011).
        A solution adaptive simulation of compressible multi-fluid flows
        with general equation of state.
        International Journal for Numerical Methods in Fluids, 67, 616-637.
    \endverbatim

SourceFiles
    JWL.C

\*---------------------------------------------------------------------------*/

#ifndef JWL_H
#define JWL_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationOfState.H"
#include "afterburnModel.H"

namespace Foam
{
namespace equationOfStates
{

/*---------------------------------------------------------------------------*\
                           Class JWL Declaration
\*---------------------------------------------------------------------------*/

class JWL
:
    public equationOfState
{
// Private data

        //- Is the phase activated
        volScalarField active_;

        //- Detonation  points
        List<vector> detonationPoints_;

        //- Detonation velocity
        dimensionedScalar vDet_;

        //- Reference pressure
        dimensionedScalar pRef_;

        //- Reference density
        dimensionedScalar rho0_;

        //- Reference Gamma
        dimensionedScalar Gamma0_;

        //- Model coefficients
        dimensionedScalar A_;
        dimensionedScalar R1_;
        dimensionedScalar B_;
        dimensionedScalar R2_;

        //- Detonation energy (initial value of e*rho)
        dimensionedScalar E0_;

        //- Afterburn model
        autoPtr<afterburnModel> afterburn_;


public:

    //- Runtime type information
    ClassName("JWL");

    // Constructor
    JWL(const volScalarField& e, const dictionary& dict);


    //- Destructor
    virtual ~JWL();


    // Member Functions

        //- Update
        virtual void update();

        //- Reference internal energy
        virtual tmp<volScalarField> e(const volScalarField& p) const;

        //- Pressure
        virtual tmp<volScalarField> pByGamma() const;

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
