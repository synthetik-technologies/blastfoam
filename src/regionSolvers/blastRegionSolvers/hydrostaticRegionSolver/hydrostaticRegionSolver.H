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
    Foam::regionSolvers::blast

Description
    Updated the hydrostatic pressure

SourceFiles
    hydrostaticRegionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef hydrostaticRegionSolver_H
#define hydrostaticRegionSolver_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionSolver.H"
#include "fluidRegionSolver.H"
#include "compressibleSystem.H"
#include "atmosphereModel.H"
#include "uniformDimensionedFields.H"

namespace Foam
{
namespace regionSolvers
{

/*---------------------------------------------------------------------------*\
                    Class hydrostaticRegionSolver Declaration
\*---------------------------------------------------------------------------*/

class hydrostatic
:
    public fluid
{
protected:
// Protected data

        //- Gravity
        uniformDimensionedVectorField g_;

        //- Compressible system
        autoPtr<compressibleSystem> fluid_;

        //- Atmosphere properties
        IOdictionary atmosphereProperties_;

        //- Atmosphere model
        autoPtr<atmosphereModel> atmosphere_;


public:

    TypeName("hydrostatic");

    // Constructors

        // Lookup compressible system type
        hydrostatic(fvMesh& mesh, const regionSolverList& regions);


    //- Destructor
    virtual ~hydrostatic();


    // Member functions

        // Edit

            //- Solve the model
            virtual void solve();

            //- Return the Courant number
            virtual scalar CoNum() const;

            //- Return the maximum Courant number
            virtual scalar maxCo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
