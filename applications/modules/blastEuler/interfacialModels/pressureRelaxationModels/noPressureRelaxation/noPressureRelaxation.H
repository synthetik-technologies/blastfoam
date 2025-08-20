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
    Foam::noPressureRelaxation

Description
    Basic ode solver for coupling phases with pressureRelaxation

SourceFiles
    noPressureRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef noPressureRelaxation_H
#define noPressureRelaxation_H

#include "pressureRelaxationSolver.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class noPressureRelaxation Declaration
\*---------------------------------------------------------------------------*/

class noPressureRelaxation
:
    public pressureRelaxationSolver
{
    // Private Member Functions

        //- Disallow copy constructor
        noPressureRelaxation(const noPressureRelaxation&);

        //- Disallow default bitwise assignment
        void operator=(const noPressureRelaxation&);

public:

    TypeName("none");

    // Constructors

        //- Construct from thermo
        noPressureRelaxation
        (
            const dictionary& dict,
            phaseSystem& fluid,
            interfacialPressureModelTable& interfacialPressureModels,
            pressureRelaxationModelTable& pressureRelaxationModels
        );


    //- Destructor
    virtual ~noPressureRelaxation();


    // Member Functions

        //- Solve the reaction system for the given time step
        virtual bool solve(const scalar& deltaT);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
