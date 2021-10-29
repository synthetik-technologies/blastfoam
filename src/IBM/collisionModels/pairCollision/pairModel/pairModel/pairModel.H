/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::pairModel

Description
    Templated pair interaction class

SourceFiles
    pairModel.C
    pairModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef immersedPairModel_H
#define immersedPairModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "pairCollisionData.H"
#include "immersedBoundaryObject.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class pairModel Declaration
\*---------------------------------------------------------------------------*/

class pairModel
{
    // Private Data

        //- The CollisionModel dictionary
        const dictionary& dict_;

        //- The coefficients dictionary
        const dictionary coeffDict_;


public:

    //- Runtime type information
    TypeName("pairModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        pairModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        pairModel(const dictionary& dict, const word& type);


    //- Destructor
    virtual ~pairModel();


    //- Selector
    static autoPtr<pairModel> New(const dictionary& dict);



    // Member Functions

        //- Return the dictionary
        const dictionary& dict() const;

        //- Return the coefficients dictionary
        const dictionary& coeffDict() const;

        //- Calculate the pair interaction between parcels
        virtual void evaluatePair
        (
            immersedBoundaryObject& OA,
            immersedBoundaryObject& OB,
            const pairCollisionData& map
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
