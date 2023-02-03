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
    Foam::displacementRelaxations::noRelaxation

Description
    No relaxation

SourceFiles
    noDisplacementRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef noDisplacementRelaxation_H
#define noDisplacementRelaxation_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "displacementRelaxation.H"

namespace Foam
{
namespace displacementRelaxations
{

/*---------------------------------------------------------------------------*\
                           Class noRelaxation Declaration
\*---------------------------------------------------------------------------*/

class noRelaxation
:
    public displacementRelaxation
{
public:

    TypeName("none");

    // Constructor
    noRelaxation(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~noRelaxation();


    // Member functions

        //- Relax the field
        virtual void relax(const label iter, pointVectorField&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace displacementRelaxations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
