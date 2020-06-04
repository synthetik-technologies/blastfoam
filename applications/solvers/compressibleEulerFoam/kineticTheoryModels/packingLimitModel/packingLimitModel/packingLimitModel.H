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
    Foam::kineticTheoryModels::packingLimitModel

SourceFiles
    packingLimitModel.C

\*---------------------------------------------------------------------------*/

#ifndef packingLimitModel_H
#define packingLimitModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "phaseModel.H"
#include "kineticTheorySystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class packingLimitModel Declaration
\*---------------------------------------------------------------------------*/

class packingLimitModel
{
    // Private data

        //- Are all particle diameters constant
        Switch constantDiameters_;


protected:

    // Protected data

        //- Reference to higher-level dictionary for re-read
        const dictionary& dict_;

        //- Const reference to kinetic theories
        const kineticTheorySystem& kt_;

        //- Const reference to mesh
        const fvMesh& mesh_;


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
            const kineticTheorySystem& kt
        ),
        (dict, kt)
    );


    // Constructors

        //- Construct from components
        packingLimitModel
        (
            const dictionary& dict,
            const kineticTheorySystem& kt
        );


    // Selectors

        static autoPtr<packingLimitModel> New
        (
            const dictionary& dict,
            const kineticTheorySystem& kt
        );


    //- Destructor
    virtual ~packingLimitModel();


    //- Disallow default bitwise copy construct
    packingLimitModel(const packingLimitModel&) = delete;
    // Member Functions

        //- Return the maximum packing limit
        tmp<volScalarField> alphaMax() const;

        // Return the packing limit for celli
        virtual scalar alphaMax
        (
            const label celli,
            const scalarList& ds
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const packingLimitModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
