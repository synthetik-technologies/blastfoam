/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::packingLimitModel

SourceFiles
    packingLimitModel.C
    packingLimitModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef packingLimitModel_H
#define packingLimitModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "phaseModel.H"
#include "masterSystem.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class packingLimitModel Declaration
\*---------------------------------------------------------------------------*/

class packingLimitModel
{
protected:

    // Protected data

        //- Const reference to kinetic theories
        const masterSystem& system_;

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Minimum max volume fraction
        scalar minAlphaMax_;


public:

    //- Runtime type information
    TypeName("packingLimitModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        packingLimitModel,
        dictionary,
        (
            const dictionary& dict,
            const masterSystem& system
        ),
        (dict, system)
    );


    // Constructors

        //- Construct from components
        packingLimitModel
        (
            const dictionary& dict,
            const masterSystem& system
        );


    // Selectors

        static autoPtr<packingLimitModel> New
        (
            const dictionary& dict,
            const masterSystem& system
        );


    //- Destructor
    virtual ~packingLimitModel();


    //- Disallow default bitwise copy construct
    packingLimitModel(const packingLimitModel&) = delete;

    // Member Functions

        //- Return the maximum packing limit
        virtual void updateAlphaMax
        (
            volScalarField::Internal& alphaMaxField
        ) const;

        // Return the packing limit for celli
        virtual scalar alphaMax
        (
            const label celli,
            SortableList<scalar>& ds,
            const bool fixedDs
        ) const = 0;

        //- Return minimum alphaMax
        inline scalar minAlphaMax() const
        {
            return minAlphaMax_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const packingLimitModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
