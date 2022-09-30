/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020 Synthetik Applied Technologies: |  Set old cell volumes after
                                refinement is done.
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
    Foam::movingAdaptiveFvMesh

SourceFiles
    movingAdaptiveFvMesh.C

Description
    Arbitary

Additional files
    adaptiveFvMesh.H
    movingAdaptiveFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef ALEFvMotionSolver_H
#define ALEFvMotionSolver_H

#include "motionSolver.H"
#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ALEFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class ALEFvMotionSolver
:
    public motionSolver,
    public fvMotionSolver
{
public:

    //- Runtime type information
    TypeName("ALE");


    // Constructors

        //- Construct from polyMesh and dictionary
        ALEFvMotionSolver
        (
            const polyMesh&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        ALEFvMotionSolver
        (
            const ALEFvMotionSolver&
        );


    //- Destructor
    ~ALEFvMotionSolver();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve()
        {}

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&)
        {}

        virtual void movePoints(const pointField&)
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ALEFvMotionSolver&) = delete;
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
