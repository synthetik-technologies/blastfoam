/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
7-2-2019 Jeff Heylmun:      Added pressureRelaxation terms
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
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
    pressureRelaxationODEI.H
    pressureRelaxationODE.C
\*---------------------------------------------------------------------------*/

#ifndef pressureRelaxationODE_H
#define pressureRelaxationODE_H

#include "ODESystem.H"
#include "ODESolver.H"
#include "volFields.H"
#include "simpleMatrix.H"
#include "phaseSystem.H"
#include "labelVector.H"
#include "pressureRelaxationModel.H"
#include "interfacialPressureModel.H"
#include "fluidBlastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class pressureRelaxationODE Declaration
\*---------------------------------------------------------------------------*/

class pressureRelaxationODE
:
    public ODESystem
{
    typedef UPtrList<phaseModel> phaseModelList;

    typedef HashTable
        <
            autoPtr<pressureRelaxationModel>,
            phasePairKey,
            phasePairKey::hash
        > pressureRelaxationModelTable;

    // Private Member Functions

        //- Disallow copy constructor
        pressureRelaxationODE(const pressureRelaxationODE&);

        //- Disallow default bitwise assignment
        void operator=(const pressureRelaxationODE&);


protected:


    // Protected data

        //- Ode dictionary
        const dictionary& dict_;

        //- Solve pressureRelaxation with ode
        Switch solvePressureRelaxation_;

        //- Reference to phase system
        const phaseSystem& fluid_;

        //- Reference to phase models
        phaseModelList phaseModels_;

        //- Local indexes of phases
        labelList phaseIndicies_;

        //- Reference to fluid thermo models
        UPtrList<fluidBlastThermo> thermos_;

        //- Reference to interfacial pressure model
        UPtrList<interfacialPressureModel> interfacialPressureModels_;

        //- PrtList of references to pressureRelaxationModels
        UPtrList<pressureRelaxationModel> pressureRelaxationModels_;

        //- Reference to phase pairs
        UPtrList<phasePair> phasePairs_;

        //- Number of equations solved
        label nEqns_;

        //- Temporary field (u1 v1 w1, u2 v2 w2 ...)
        scalarField q_;

        //- Temporary rate-of-change of variables field
        scalarField dqdt_;

        //- Stored time step
        volScalarField::Internal deltaT_;

        //- Ode solver
        autoPtr<ODESolver> odeSolver_;


public:


    // Constructors

        //- Construct from thermo
        pressureRelaxationODE(phaseSystem& fluid, pressureRelaxationModelTable& pressureRelaxationModels);


    //- Destructor
    virtual ~pressureRelaxationODE();


    // Member Functions

        bool solve() const
        {
            return solvePressureRelaxation_;
        }

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
