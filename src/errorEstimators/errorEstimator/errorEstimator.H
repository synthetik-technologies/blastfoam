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
    Foam::errorEstimator

Description
    Base class used to estimate error within a cell/across faces

SourceFiles
    errorEstimator.C
    newErrorEstimator.C

\*---------------------------------------------------------------------------*/

#ifndef errorEstimator_H
#define errorEstimator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class errorEstimator Declaration
\*---------------------------------------------------------------------------*/

class errorEstimator
:
    public volScalarField
{
protected:
    // Protected data

        //- Reference to mesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    ClassName("errorEstimator");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            errorEstimator,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );

    // Constructor
    errorEstimator(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~errorEstimator();


    // Selectors

        static autoPtr<errorEstimator> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        //- Update
        virtual void update() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
