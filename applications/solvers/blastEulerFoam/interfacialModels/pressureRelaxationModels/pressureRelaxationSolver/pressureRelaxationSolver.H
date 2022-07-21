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
    Foam::pressureRelaxationSolver

Description
    Base class for pressure relaxation

SourceFiles
    pressureRelaxationSolver.C

\*---------------------------------------------------------------------------*/

#ifndef pressureRelaxationSolver_H
#define pressureRelaxationSolver_H

#include "volFields.H"
#include "simpleMatrix.H"
#include "phaseSystem.H"
#include "labelVector.H"
#include "pressureRelaxationModel.H"
#include "interfacialPressureModel.H"
#include "fluidBlastThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class pressureRelaxationSolver Declaration
\*---------------------------------------------------------------------------*/

class pressureRelaxationSolver
{
public:

    typedef UPtrList<phaseModel> phaseModelList;

    typedef HashTable
        <
            autoPtr<interfacialPressureModel>,
            phasePairKey,
            phasePairKey::hash
        > interfacialPressureModelTable;

    typedef HashTable
        <
            autoPtr<pressureRelaxationModel>,
            phasePairKey,
            phasePairKey::hash
        > pressureRelaxationModelTable;

private:
    // Private Member Functions

        //- Disallow copy constructor
        pressureRelaxationSolver(const pressureRelaxationSolver&);

        //- Disallow default bitwise assignment
        void operator=(const pressureRelaxationSolver&);


protected:


    // Protected data

        //- Solve pressureRelaxation with ode
        Switch solvePressureRelaxation_;

        //- Reference to phase system
        const phaseSystem& fluid_;

        //- Reference to phase models
        mutable phaseModelList phaseModels_;

        //- Local indexes of phases
        Map<label> phaseIndicies_;

        //- Reference to fluid thermo models
        mutable UPtrList<fluidBlastThermo> thermos_;

        //- Reference to interfacial pressure model
        UPtrList<interfacialPressureModel> interfacialPressureModels_;

        //- PrtList of references to pressureRelaxationModels
        UPtrList<pressureRelaxationModel> pressureRelaxationModels_;

        //- Number of equations solved
        label nEqns_;

public:

    //- Runtime type information
    TypeName("pressureRelaxationSolver");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            pressureRelaxationSolver,
            dictionary,
            (
                phaseSystem& fluid,
                interfacialPressureModelTable& interfacialPressureModels,
                pressureRelaxationModelTable& pressureRelaxationModels
            ),
            (fluid, interfacialPressureModels, pressureRelaxationModels)
        );

    // Constructors

        //- Construct from thermo
        pressureRelaxationSolver
        (
            phaseSystem& fluid,
            interfacialPressureModelTable& interfacialPressureModels,
            pressureRelaxationModelTable& pressureRelaxationModels,
            const bool needIPModels,
            const bool needPRModels
        );

        //- Construct from thermo
        pressureRelaxationSolver
        (
            phaseSystem& fluid,
            interfacialPressureModelTable& interfacialPressureModels,
            pressureRelaxationModelTable& pressureRelaxationModels
        );


    //- Destructor
    virtual ~pressureRelaxationSolver();

    // Selectors

        static autoPtr<pressureRelaxationSolver> New
        (
            phaseSystem& fluid,
            interfacialPressureModelTable& interfacialPressureModels,
            pressureRelaxationModelTable& pressureRelaxationModels
        );

    // Member Functions

        bool solve() const
        {
            return solvePressureRelaxation_;
        }

        //- Solve the reaction system for the given time step
        //  and return the characteristic time
        virtual scalar solve(const scalar& deltaT)
        {
            return deltaT;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
