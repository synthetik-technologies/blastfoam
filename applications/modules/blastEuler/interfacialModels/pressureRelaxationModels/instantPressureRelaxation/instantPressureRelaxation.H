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
#include "univariateRootSolver.H"
#include "IntegratorsFwd.H"
#include "EquationsFwd.H"
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

        class alphaPIIntegrator
        :
            public ScalarEquation
        {
            mutable phaseModelList phases_;
            mutable UPtrList<fluidBlastThermo> thermos_;

        public:
            alphaPIIntegrator
            (
                phaseModelList& phases,
                UPtrList<fluidBlastThermo>& thermos
            )
            :
                ScalarEquation(0.0, 1.0),
                phases_(phases),
                thermos_(thermos)
            {}

            virtual ~alphaPIIntegrator() {}

            //- Return the number of derivatives
            virtual label nDerivatives() const
            {
                return 0;
            }

            virtual scalar fx
            (
                const scalar x,
                const label celli
            ) const
            {
                // const fluidBlastThermo& thermo1 = thermos_[0];
                // const fluidBlastThermo& thermo2 = thermos_[1];
                //
                // const scalar alpha1 = max(min(1.0, x), 0.0);
                // const scalar alpha2 = 1.0 - alpha1;
                //
                // return
                //     alpha1*thermo1.p()[celli]
                //   + alpha2*thermo2.p()[celli];

                phaseModel& phase1 = phases_[0];
                phaseModel& phase2 = phases_[1];
                fluidBlastThermo& thermo1 = thermos_[0];
                fluidBlastThermo& thermo2 = thermos_[1];

                const scalar alpha1 = x;
                const scalar alpha2 = 1.0 - alpha1;

                const scalar alphaRho1 = phase1.alphaRho()[celli];
                const scalar alphaRho2 = phase2.alphaRho()[celli];
                if (alphaRho1 < 1e-6)
                {
                    thermo2.rhoRef()[celli] = alphaRho2/max(alpha2, 1e-6);
                    return thermo2.cellpRhoT(celli);
                }
                else if (alphaRho2 < 1e-6)
                {
                    thermo1.rhoRef()[celli] = alphaRho1/max(alpha1, 1e-6);
                    return thermo1.cellpRhoT(celli);
                }

                thermo1.rhoRef()[celli] = alphaRho1/max(alpha1, 1e-6);
                thermo2.rhoRef()[celli] = alphaRho2/max(alpha2, 1e-6);
                return
                    alpha1*thermo1.cellpRhoT(celli)
                  + alpha2*thermo2.cellpRhoT(celli);


            }
        };

        class simplePIEqn
        :
            public ScalarEquation
        {
            mutable phaseModelList phases_;
            mutable UPtrList<fluidBlastThermo> thermos_;
            const scalarField& alpha0_;
            const scalarField& e0_;
            const scalarField& Pi0_;

        public:
            simplePIEqn
            (
                phaseModelList& phases,
                UPtrList<fluidBlastThermo>& thermos,
                const scalarField& alpha0,
                const scalarField& e0,
                const scalarField& Pi0

            )
            :
                ScalarEquation(0, 1),
                phases_(phases),
                thermos_(thermos),
                alpha0_(alpha0),
                e0_(e0),
                Pi0_(Pi0)
            {}

            virtual ~simplePIEqn() {}

            //- Return the number of derivatives
            virtual label nDerivatives() const
            {
                return 1;
            }

            virtual scalar fx
            (
                const scalar x,
                const label celli
            ) const
            {
                phaseModel& phase1 = phases_[0];
                phaseModel& phase2 = phases_[1];
                fluidBlastThermo& thermo1 = thermos_[0];
                fluidBlastThermo& thermo2 = thermos_[1];

                const scalar alpha1 = x;
                const scalar alpha2 = 1.0 - alpha1;

                // if (alpha1 < 1e-6)
                // {
                //     return 0.0;
                // }
                // else if (alpha2 < 1e-6)
                // {
                //     return 0.0;
                // }

                const scalar alphaRho1 = phase1.alphaRho()[celli];
                const scalar alphaRho2 = phase2.alphaRho()[celli];

                const scalar PI1 = 0.5*(thermo1.p()[celli] + Pi0_[0]);
                const scalar PI2 = 0.5*(thermo2.p()[celli] + Pi0_[1]);
                const scalar PI = alpha1*PI1 + alpha2*PI2;


                thermo1.rhoRef()[celli] = alphaRho1/max(alpha1, 1e-6);
                thermo1.he()[celli] =
                    e0_[0] - PI*(alpha1 - alpha0_[0])/max(alphaRho1, 1e-6);
                thermo1.T()[celli] = thermo1.cellThe
                (
                    thermo1.he()[celli],
                    thermo1.T()[celli],
                    celli
                );
                thermo1.p()[celli] = thermo1.cellpRhoT(celli);

                thermo2.rhoRef()[celli] = alphaRho2/max(alpha2, 1e-6);
                thermo2.he()[celli] =
                    e0_[1] + PI*(alpha1 - alpha0_[0])/max(alphaRho1, 1e-6);
                thermo2.T()[celli] = thermo2.cellThe
                (
                    thermo2.he()[celli],
                    thermo2.T()[celli],
                    celli
                );
                thermo2.p()[celli] = thermo2.cellpRhoT(celli);

                return thermo1.cellpRhoT(celli) - thermo2.cellpRhoT(celli);
            }

            virtual scalar dfdx
            (
                const scalar x,
                const label celli
            ) const
            {
                phaseModel& phase1 = phases_[0];
                phaseModel& phase2 = phases_[1];
                fluidBlastThermo& thermo1 = thermos_[0];
                fluidBlastThermo& thermo2 = thermos_[1];

                const scalar alpha1 = max(min(1.0, x), 0.0);
                const scalar alpha2 = 1.0 - alpha1;

                const scalar alphaRho1 = phase1.alphaRho()[celli];
                const scalar alphaRho2 = phase2.alphaRho()[celli];

                const scalar PI =
                    alpha1*thermo1.p()[celli] + alpha2*thermo2.p()[celli];

                thermo1.rhoRef()[celli] = alphaRho1/max(alpha1, 1e-6);
                thermo1.he()[celli] =
                    e0_[0] - PI*(alpha1 - alpha0_[0])/max(alphaRho1, 1e-10);
                thermo1.T()[celli] = thermo1.cellThe
                (
                    thermo1.he()[celli],
                    thermo1.T()[celli],
                    celli
                );

                thermo2.rhoRef()[celli] = alphaRho2/max(alpha2, 1e-6);
                thermo2.he()[celli] =
                    e0_[1] - PI*(alpha2 - alpha0_[1])/max(alphaRho2, 1e-10);
                thermo2.T()[celli] = thermo2.cellThe
                (
                    thermo2.he()[celli],
                    thermo2.T()[celli],
                    celli
                );

                thermo1.p()[celli] = thermo1.cellpRhoT(celli);
                thermo2.p()[celli] = thermo2.cellpRhoT(celli);


                const scalar dp1dalpha =
                  - thermo1.celldpde(celli)*PI/max(alphaRho1, 1e-10)
                  - thermo1.celldpdRho(celli)*alphaRho1/sqr(max(alpha1, 1e-6));
                const scalar dp2dalpha =
                    thermo2.celldpde(celli)*PI/max(alphaRho2, 1e-10)
                  + thermo2.celldpdRho(celli)*alphaRho2/sqr(max(alpha2, 1e-6));

                return dp1dalpha - dp2dalpha;
            }
        };
protected:


    // Protected data

        //- alphaPI integration equation
        alphaPIIntegrator alphaPIEqn_;

        //- alphaPI integration equation
        simplePIEqn simpleEqn_;

        //- Root solver
        autoPtr<rootSolver> rootSolver_;

        //- Root solver
        autoPtr<univariateRootSolver> uniRootSolver_;

        //- Root solver
        autoPtr<scalarIntegrator> integrator_;

        //- Saved old values
        scalarField alpha0_;
        scalarField e0_;
        scalarField Pi0_;

        //- Procedure 1 relaxation
        void relax1(const scalar deltaT);

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
                const ScalarMultivariateEquation::VarType& x,
                const label li,
                scalarList& fx
            ) const;

            virtual void jacobian
            (
                const ScalarMultivariateEquation::VarType& x,
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
