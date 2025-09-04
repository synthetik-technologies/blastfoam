/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022-2025
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
    Foam::instantPressureRelaxation

Description
    Using a multivariate root solver find the volume fractions and energy
    that are in equilibrium between all phases

SourceFiles
    instantPressureRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef instantPressureRelaxation_H
#define instantPressureRelaxation_H

#include "pressureRelaxationSolver.H"
#include "rootSolver.H"
#include "MultivariateEquationsFwd.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class pressureRelaxationODE Declaration
\*---------------------------------------------------------------------------*/

class instantPressureRelaxation
:
    public pressureRelaxationSolver,
    public ScalarMultivariateEquation
{
    // Private Member Functions

        //- Disallow copy constructor
        instantPressureRelaxation(const instantPressureRelaxation&);

        //- Disallow default bitwise assignment
        void operator=(const instantPressureRelaxation&);

protected:


    // Protected data

        //- Root solver
        autoPtr<rootSolver> rootSolver_;

        //- Saved old values
        DynamicList<scalar> alpha0_;
        DynamicList<scalar> alphaRho0_;
        DynamicList<scalar> rho0_;
        DynamicList<scalar> e0_;
        scalar PI0_;
        scalar alphaMax_;

        DynamicList<label> phases_;

        //- Use multivariate root solver
        bool useRootSolver_;

        //- Procedure 4 relaxation
        void relax4(const scalar deltaT);


public:

    TypeName("instant");

    //- Constructor
        instantPressureRelaxation
        (
            const dictionary& dict,
            phaseSystem& fluid,
            interfacialPressureModelTable& interfacialPressureModels,
            pressureRelaxationModelTable& pressureRelaxationModels
        );


    //- Destructor
    virtual ~instantPressureRelaxation();


    // Member Functions

        //- Return if instant relaxation
        virtual bool instant() const
        {
            return true;
        }

        //- Solve the reaction system for the given time step
        virtual bool solve(const scalar& deltaT);


    // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            virtual label nVars() const
            {
                return pressureRelaxationSolver::nEqns_;
            }
            virtual label nEqns() const
            {
                return pressureRelaxationSolver::nEqns_;
            }

            //- Number of implemented derivatives
            virtual label nDerivatives() const
            {
                return 0;
            }

            virtual void FX
            (
                const ScalarMultivariateEquation::VarType& x,
                const label li,
                scalarList& fx
            ) const;

        // //- Calculate the Jacobian of the equation
        // virtual void jacobian
        // (
        //     const typename ScalarMultivariateEquation::VarType& x,
        //     const label li,
        //     scalarList& fx,
        //     RectangularMatrix<scalar>& J
        // ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
