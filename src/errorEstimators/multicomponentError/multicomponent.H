/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2021
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
    Foam::errorEstimators::multicomponent

Description
    Error is defined by multiple errorEstimators

SourceFiles
    multicomponent.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponent_H
#define multicomponent_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "errorEstimator.H"

namespace Foam
{
namespace errorEstimators
{

/*---------------------------------------------------------------------------*\
                           Class multicomponent Declaration
\*---------------------------------------------------------------------------*/

class multicomponent
:
    public errorEstimator
{

    //- Names of error estimators
    wordList names_;

    //- List of error estimators
    PtrList<errorEstimator> errors_;

    //- Maximum refinement level for the estimator
    mutable labelList maxRefinement_;


public:

    //- Runtime type information
    TypeName("multicomponent");

    // Constructor
    multicomponent
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& name
    );


    //- Destructor
    virtual ~multicomponent();


    // Member Functions

        //- Update error
        virtual void update(const bool scale = true);

        //- Return non constant reference to error field
        virtual labelList maxRefinement() const;

        //- Read level values
        virtual void read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace errorEstimators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
