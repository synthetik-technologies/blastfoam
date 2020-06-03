/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
    Foam::activationModels::reactionRateActivation

Description
    Constant activation rate

SourceFiles
    reactionRateActivation.C

\*---------------------------------------------------------------------------*/

#ifndef reactionRateActivation_H
#define reactionRateActivation_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "activationModel.H"

namespace Foam
{
namespace activationModels
{

/*---------------------------------------------------------------------------*\
                           Class reactionRateActivation Declaration
\*---------------------------------------------------------------------------*/

class reactionRateActivation
:
    public activationModel
{
// Private data

    //- Activation rate
    dimensionedScalar w_;

    //- Return the time rate of chage of lambda
    virtual tmp<volScalarField> delta() const;


public:

    //- Runtime type information
    TypeName("reactionRate");

    // Constructor
    reactionRateActivation
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );


    //- Destructor
    virtual ~reactionRateActivation();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace activationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
