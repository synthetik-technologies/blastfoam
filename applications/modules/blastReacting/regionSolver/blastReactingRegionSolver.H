/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::regionSolvers::blastReacting

Description
    Solver for a reacting fluid phase using standard OpenFOAM thermodynamics

SourceFiles
    blastReactingRegionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef blastReactingRegionSolver_H
#define blastReactingRegionSolver_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "blastRegionSolver.H"
#include "reactingCompressibleSystem.H"
#include "uniformDimensionedFields.H"

namespace Foam
{
namespace regionSolvers
{

/*---------------------------------------------------------------------------*\
                    Class reactingRegionSolver Declaration
\*---------------------------------------------------------------------------*/

class blastReacting
:
    public blast
{
public:

    TypeName("blastReacting");

    // Constructor
    blastReacting
    (
        fvMesh& mesh,
        const regionSolverList& regions
    );

    //- Destructor
    virtual ~blastReacting();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
