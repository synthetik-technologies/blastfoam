/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
02-06-2020  Jeff Heylmun    : Modified ODE system to solve radiation
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
    Foam::radiationODE

Description
    Basic ode solver for relaxing temperatures using radiation heat transfer.

SourceFiles
    radiationODEI.H
    radiationODE.C
\*---------------------------------------------------------------------------*/

#ifndef radiationODE_H
#define radiationODE_H

#include "ODESystem.H"
#include "ODESolver.H"
#include "volFields.H"
#include "simpleMatrix.H"
#include "basicThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class radiationModel;

/*---------------------------------------------------------------------------*\
                      Class radiationODE Declaration
\*---------------------------------------------------------------------------*/

class radiationODE
:
    public ODESystem
{
    // Private data

        //- Reference to radiation model
        const radiationModel& rad_;

        //- Reference to thermodynamic model
        const basicThermoModel& thermo_;


    // Private Member Functions


        //- Disallow copy constructor
        radiationODE(const radiationODE&);

        //- Disallow default bitwise assignment
        void operator=(const radiationODE&);


protected:


    // Protected data

        //- Ode dictionary
        const dictionary& dict_;

        //- Is the ode solved
        bool solve_;

        //- Number of equations solved
        label nEqns_;

        //- Temporary field (u1 v1 w1, u2 v2 w2 ...)
        scalarField q_;

        //- Temporary rate-of-change of variables field
        scalarField dqdt_;

        //- Stored time step
        scalarField deltaT_;

        //- Current cell
        label celli_ = 0;

        //- Ode solver
        autoPtr<ODESolver> odeSolver_;


public:


    // Constructors

        //- Construct from thermo
        radiationODE
        (
            const radiationModel& rad,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~radiationODE();


    // Member Functions

        bool solve() const
        {
            return solve_;
        }

        //- Solve the reaction system for the given time step
        //  and return the characteristic time
        virtual scalar solve
        (
            const scalar& deltaT,
            volScalarField& rhoE
        );

    // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            virtual label nEqns() const
            {
                return nEqns_;
            }

            virtual void derivatives
            (
                const scalar t,
                const scalarField& q,
                scalarField& dqdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& q,
                scalarField& dqdt,
                scalarSquareMatrix& J
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
