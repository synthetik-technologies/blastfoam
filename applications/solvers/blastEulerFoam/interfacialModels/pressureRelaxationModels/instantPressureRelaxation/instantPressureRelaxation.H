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

        //- Root solver dictionary
        const dictionary& dict_;

        //- Root solver
        autoPtr<rootSolver> rootSolver_;

        //- Original densitites
        scalarList rho0_;

        //- Original energies
        scalarList e0_;

        //- Original interface pressure
        scalar PI0_;


public:

    TypeName("instant");

    //- Constructor
        instantPressureRelaxation
        (
            phaseSystem& fluid,
            interfacialPressureModelTable& interfacialPressureModels,
            pressureRelaxationModelTable& pressureRelaxationModels
        );


    //- Destructor
    virtual ~instantPressureRelaxation();


    // Member Functions

        //- Solve the reaction system for the given time step
        //  and return the characteristic time
        virtual scalar solve(const scalar& deltaT);


    // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            virtual label nEqns() const
            {
                return pressureRelaxationSolver::nEqns_;
            }

            //- Number of implemented derivatives
            virtual label nDerivatives() const
            {
                return 1;
            }

            virtual void FX
            (
                const scalarList& x,
                const label li,
                scalarList& fx
            ) const;

            virtual void jacobian
            (
                const scalarList& x,
                const label li,
                scalarList& fx,
                RectangularMatrix<scalar>& J
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
