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
    Foam::MGEquationOfState

Description

    Simple class to transfer the correct quantities for MieGruniesen based
    equation of state

SourceFiles
    MGEquationOfState.C

\*---------------------------------------------------------------------------*/

#ifndef MGEquationOfState_H
#define MGEquationOfState_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"
#include "CpCvModel.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class MGEquationOfState Declaration
\*---------------------------------------------------------------------------*/

class MGEquationOfState
{
protected:

    //- Pointer to cpcv model
    const CpCvModel* cpcvPtr_;

    //- Return the MieGruniesen coefficient
    virtual scalar GammaMG
    (
        const scalar rho,
        const scalar e
    ) const = 0;

    //- Return the deviation from ideal pressure
    virtual scalar Pi
    (
        const scalar rho,
        const scalar e
    ) const = 0;

    //- Return the nonideal speed of sound contribution
    virtual scalar delta
    (
        const scalar p,
        const scalar rho,
        const scalar e
    ) const = 0;


public:

        virtual ~MGEquationOfState()
        {}

    // Member Functions

        //- Is temperature used
        bool temperatureBased() const
        {
            return false;
        }

        //- Set the pointer to the CpCvModel
        void set(const CpCvModel& cpcv)
        {
            cpcvPtr_ = &cpcv;
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
            return
            limit
          ? max((GammaMG(rho, e) - 1.0)*rho*e - Pi(rho, e), 0.0)
          : (GammaMG(rho, e) - 1.0)*rho*e - Pi(rho, e);
        }

        //- Return the speed of sound
        virtual scalar cSqr
        (
            const scalar p,
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar Cv
        ) const
        {
            scalar G(GammaMG(rho, e));
            scalar h((G*p + Pi(rho, e))/((G - 1.0)*max(rho, 1e-10)));
            return (h - delta(p, rho, e))*(G - 1.0);
        }

        //- Return derivative of pressure w.r.t. Temperature
        inline scalar dpdT
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const
        {
            NotImplemented;
            return 0.0;
        }

        //- Return the MieGruniesen coefficient
        inline scalar Gamma
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar Cv
        ) const
        {
            return this->GammaMG(rho, e);
        }

        //- Heat capacity at constant volume [J/kg/K]
        virtual scalar Cv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const
        {
            return 0.0;
        }

        //- Heat capacity at constant pressure [J/kg/K]
        virtual scalar Cp
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const
        {
            return 0.0;
        }

        //- Heat capacity difference [J/kg/K]
        virtual scalar CpMCv
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return (this->GammaMG(rho, e) - 1.0)*cpcvPtr_->Cv(rho, e, T);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
