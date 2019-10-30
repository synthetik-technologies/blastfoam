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
    Foam::multiphaseEquationOfState

Description
    Class to calculate mixture properties of a collection of more than two
    equation of states.

    References:
    \verbatim
        Zheng, H.W., Shu, C., Chew, Y.T., Qin, N.  (2011).
        A solution adaptive simulation of compressible multi-fluid flows
        with general equation of state.
        International Journal for Numerical Methods in Fluids, 67, 616-637.
    \endverbatim

SourceFiles
    multiphaseEquationOfState.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseEquationOfState_H
#define multiphaseEquationOfState_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationOfState.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multiphaseEquationOfState Declaration
\*---------------------------------------------------------------------------*/

class multiphaseEquationOfState
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

        //- Equation of state for phases
        PtrList<equationOfState> eos_;

        //- Volume fraction of phases
        UPtrList<volScalarField> alphas_;

        //- Volume weighted density of phases
        UPtrList<volScalarField> rhos_;

        //- Volume weighted density of phases
        UPtrList<volScalarField> alphaRhos_;

        //- Volume fraction fluxes
        UPtrList<surfaceScalarField> alphaPhis_;

        //- Mass fluxes
        UPtrList<surfaceScalarField> alphaRhoPhis_;


public:


    // Constructor
    multiphaseEquationOfState
    (
        volScalarField& rho,
        volScalarField& e,
        volScalarField& p,
        const dictionary& dict
    );


    //- Destructor
    virtual ~multiphaseEquationOfState();


    // Member Functions

        //- Update pressure
        void updateP();

        //- Return shared speed of sound
        tmp<volScalarField> c() const;


    // Access functions

        //- Const reference to volume fraction
        const UPtrList<volScalarField> alphas() const
        {
            return alphas_;
        }

        //- Non-const reference to volume fraction
        UPtrList<volScalarField>& alphas()
        {
            return alphas_;
        }

        //- Non-const reference to phase densities
        UPtrList<volScalarField>& rhos()
        {
            return rhos_;
        }

        //- Const reference to volume weighted density of phases
        const UPtrList<volScalarField> alphaRhos() const
        {
            return alphaRhos_;
        }

        //- Non-const reference to volume weighted density of phases
        UPtrList<volScalarField>& alphaRhos()
        {
            return alphaRhos_;
        }

        //- Const reference to phase densities
        const UPtrList<surfaceScalarField> alphaPhis() const
        {
            return alphaPhis_;
        }

        //- Non-const reference to phase densities
        UPtrList<surfaceScalarField>& alphaPhis()
        {
            return alphaPhis_;
        }

        //- Const reference to volume weighted density of phases
        const UPtrList<surfaceScalarField> alphaRhoPhis() const
        {
            return alphaRhoPhis_;
        }

        //- Non-const reference to volume weighted density of phases
        UPtrList<surfaceScalarField>& alphaRhoPhis()
        {
            return alphaRhoPhis_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
