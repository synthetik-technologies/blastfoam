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
    Foam::twoPhaseEquationOfState

Description
    Class to calculate mixture properties of two equation of states

    References:
    \verbatim
        Zheng, H.W., Shu, C., Chew, Y.T., Qin, N.  (2011).
        A solution adaptive simulation of compressible multi-fluid flows
        with general equation of state.
        International Journal for Numerical Methods in Fluids, 67, 616-637.
    \endverbatim

SourceFiles
    twoPhaseEquationOfState.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseEquationOfState_H
#define twoPhaseEquationOfState_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationOfState.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class twoPhaseEquationOfState Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseEquationOfState
{
// Private data

        //- Names of phases
        wordList phases_;

        //- Density field
        volScalarField& rho_;

        //- Internal energy
        volScalarField& e_;

        //- Shared pressure
        volScalarField& p_;

        //- Equation of state for phase 1
        autoPtr<equationOfState> eos1_;

        //- Volume fraction
        volScalarField& alpha_;

        //- Volume weighted density of phase 1
        volScalarField& alphaRho1_;

        //- Equation of state for phase 2
        autoPtr<equationOfState> eos2_;

        //- Volume weighted density of phase 2
        volScalarField& alphaRho2_;


public:


    // Constructor
    twoPhaseEquationOfState
    (
        volScalarField& rho,
        volScalarField& e,
        volScalarField& p,
        const dictionary& dict
    );


    //- Destructor
    virtual ~twoPhaseEquationOfState();


    // Member Functions

        //- Update pressure
        void updateP();

        //- Return shared speed of sound
        tmp<volScalarField> c() const;

        //- Const reference to volume fraction
        const volScalarField& alpha() const
        {
            return alpha_;
        }

        //- Non-const reference to volume fraction
        volScalarField& alpha()
        {
            return alpha_;
        }

        //- Const reference to phase 1 density
        const volScalarField& rho1() const
        {
            return eos1_->rho();
        }

        //- Non-const reference to phase 1 density
        volScalarField& rho1()
        {
            return eos1_->rho();
        }

        //- Const reference to volume weighted density of phase 1
        const volScalarField& alphaRho1() const
        {
            return alphaRho1_;
        }

        //- Non-const reference to volume weighted density of phase 1
        volScalarField& alphaRho1()
        {
            return alphaRho1_;
        }

        //- Const reference to volume weighted density of phase 1
        const surfaceScalarField& alphaPhi() const
        {
            return eos1_->alphaPhi();
        }

        //- Non-const reference to volume weighted density of phase 1
        surfaceScalarField& alphaPhi()
        {
            return eos1_->alphaPhi();
        }

        //- Const reference to volume weighted density of phase 1
        const surfaceScalarField& alphaRhoPhi1() const
        {
            return eos1_->alphaRhoPhi();
        }

        //- Non-const reference to volume weighted density of phase 1
        surfaceScalarField& alphaRhoPhi1()
        {
            return eos1_->alphaRhoPhi();
        }

        //- Const reference to phase 2 density
        const volScalarField& rho2() const
        {
            return eos2_->rho();
        }

        //- Non-const reference to phase 2 density
        volScalarField& rho2()
        {
            return eos2_->rho();
        }

        //- Const reference to volume weighted density of phase 2
        const volScalarField& alphaRho2() const
        {
            return alphaRho2_;
        }

        //- Non-const reference to volume weighted density of phase 2
        volScalarField& alphaRho2()
        {
            return alphaRho2_;
        }

        //- Const reference to volume weighted density of phase 2
        const surfaceScalarField& alphaRhoPhi2() const
        {
            return eos2_->alphaRhoPhi();
        }

        //- Non-const reference to volume weighted density of phase 2
        surfaceScalarField& alphaRhoPhi2()
        {
            return eos2_->alphaRhoPhi();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
