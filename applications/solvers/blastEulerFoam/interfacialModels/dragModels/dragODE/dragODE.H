/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
7-2-2019 Jeff Heylmun:      Added drag terms
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
    Foam::dragODE
Description
    Basic ode solver for coupling two phases with drag
SourceFiles
    dragODEI.H
    dragODE.C
\*---------------------------------------------------------------------------*/

#ifndef dragODE_H
#define dragODE_H

#include "ODESystem.H"
#include "ODESolver.H"
#include "volFields.H"
#include "simpleMatrix.H"
#include "phaseSystem.H"
#include "labelVector.H"
#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class dragODE Declaration
\*---------------------------------------------------------------------------*/

class dragODE
:
    public ODESystem
{
    typedef PtrListDictionary<phaseModel> phaseModelList;
    typedef HashTable
        <
            autoPtr<phasePair>,
            phasePairKey,
            phasePairKey::hash
        > phasePairTable;
    typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<dragModel>>,
            phasePairKey,
            phasePairKey::hash
        > dragModelTable;

    // Private Member Functions

        //- Update velocities in phase models
        void setUs();

        //- Update q vector
        void setq();

        //- Disallow copy constructor
        dragODE(const dragODE&);

        //- Disallow default bitwise assignment
        void operator=(const dragODE&);


protected:


    // Protected data

        //- Ode dictionary
        const dictionary& dict_;

        //- Solve drag with ode
        Switch solveDrag_;

        //- Reference to phase system
        const phaseSystem& fluid_;

        //- Reference to phase models
        phaseModelList& phaseModels_;

        //- PrtList of references to dragModels
        UPtrList<BlendedInterfacialModel<dragModel>> dragModels_;

        //- Reference to phase pairs
        UPtrList<phasePair> phasePairs_;

        //- vector of solution dimensions
        labelVector solutionD_;

        //- number of solution dimensions
        label nDims_;

        //- Starting index of the phases
        labelList startI_;

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

        //- Calculate the number of equations to solve
        //  and set the starting index of the phases
        label calcNEqns();


public:


    // Constructors

        //- Construct from thermo
        dragODE(phaseSystem& fluid, dragModelTable& dragModels);


    //- Destructor
    virtual ~dragODE();


    // Member Functions

        bool solveDrag() const
        {
            return solveDrag_;
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
