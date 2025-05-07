/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::pressureRelaxationODE

Description
    Basic ode solver for coupling phases with pressureRelaxation

SourceFiles
    pressureRelaxationODE.C

\*---------------------------------------------------------------------------*/

#ifndef pressureRelaxationODE_H
#define pressureRelaxationODE_H

#include "pressureRelaxationSolver.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "volFields.H"
#include "simpleMatrix.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class pressureRelaxationODE Declaration
\*---------------------------------------------------------------------------*/

class pressureRelaxationODE
:
    public pressureRelaxationSolver,
    public ODESystem
{
    // Private Member Functions

        //- Disallow copy constructor
        pressureRelaxationODE(const pressureRelaxationODE&);

        //- Disallow default bitwise assignment
        void operator=(const pressureRelaxationODE&);


protected:
    // Protected data

        //- Ode dictionary
        const dictionary& dict_;

        //- Temporary field (u1 v1 w1, u2 v2 w2 ...)
        scalarField q_;

        //- Temporary rate-of-change of variables field
        scalarField dqdt_;

        //- Stored time step
        volScalarField::Internal deltaT_;

        //- Ode solver
        autoPtr<ODESolver> odeSolver_;


public:

    TypeName("ode");

    // Constructors

        //- Construct from thermo
        pressureRelaxationODE
        (
            phaseSystem& fluid,
            interfacialPressureModelTable& interfacialPressureModels,
            pressureRelaxationModelTable& pressureRelaxationModels
        );


    //- Destructor
    virtual ~pressureRelaxationODE();


    // Member Functions

        //- Solve the reaction system for the given time step
        //  and return the characteristic time
        virtual scalar solve(const scalar& deltaT);

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
                const label li,
                scalarField& dqdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& q,
                const label li,
                scalarField& dqdt,
                scalarSquareMatrix& J
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
