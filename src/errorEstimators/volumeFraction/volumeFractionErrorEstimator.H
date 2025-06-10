/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023
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
    Foam::errorEstimators::volumeFraction

Description
    Error is the value of the given field

SourceFiles
    volumeFractionErrorEstimator.C

\*---------------------------------------------------------------------------*/

#ifndef volumeFractionErrorEstimator_H
#define volumeFractionErrorEstimator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fieldValueErrorEstimator.H"

namespace Foam
{
namespace errorEstimators
{

/*---------------------------------------------------------------------------*\
                           Class volumeFraction Declaration
\*---------------------------------------------------------------------------*/

class volumeFraction
:
    public fieldValue
{
public:

    //- Runtime type information
    TypeName("volumeFraction");

    // Constructor
    volumeFraction(const fvMesh& mesh, const dictionary& dict, const word& name);


    //- Destructor
    virtual ~volumeFraction();


    // Member Functions

        //- Read level values
        virtual void read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace errorEstimators
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
