/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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
    Foam::equationOfState

Description
    Simple class to transfer the correct quantities for a temperature based
    equation of state

SourceFiles
    equationOfState.C

\*---------------------------------------------------------------------------*/

#ifndef equationOfState_H
#define equationOfState_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UautoPtr.H"
#include "CpCvModel.H"

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class equationOfState Declaration
\*---------------------------------------------------------------------------*/

class equationOfState
{
protected:

    virtual scalar pRhoT(const scalar rho, const scalar T) const = 0;


public:

    equationOfState()
    {}

    virtual ~equationOfState()
    {}


    // Member Functions

        //- Is temperature used
        bool temperatureBased() const
        {
            return true;
        }

        //- Return pressure
        virtual scalar p
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const bool limit = true
        ) const
        {
            return limit ? max(this->pRhoT(rho, T), 0.0) : this->pRhoT(rho, T);
        }

        //- Return the speed of sound
        virtual scalar cSqr
        (
            const scalar p,
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
        ) const
        {
            return
                (sqr(dpdT(rho, e, T))*T/cv - dpdv(rho, e, T))
               /sqr(max(rho, 1e-10));
        }

        //- Return the MieGruniesen coefficient
        virtual scalar Gamma
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
        ) const
        {
            return max(dpdT(rho, e, T)/(cv*max(rho, 1e-10)), small) + 1.0;
        }

        //- Return derivative of pressure w.r.t. specific volume
        virtual scalar dpdv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const = 0;

        //- Return derivative of pressure w.r.t. Temperature
        virtual scalar dpdT
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const = 0;

        //- Return derivative of pressure w.r.t. internal energy
        scalar dpde
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const
        {
            NotImplemented;
            return 0.0;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
