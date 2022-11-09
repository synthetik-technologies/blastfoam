/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::rainfallModel

Description
    Generic class for finding a root of a univariate equation

SourceFiles
    rainfallModel.C

\*---------------------------------------------------------------------------*/

#ifndef rainfallModel_H
#define rainfallModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class rainfallModel Declaration
\*---------------------------------------------------------------------------*/

class rainfallModel
{
protected:

    //- Reference to the mesh
    const fvMesh& mesh_;

    //- Start time for rain
    scalar start_;

    //- End time for rain
    scalar end_;


public:

    //- Runtime type information (debug defined in rootSolver)
    TypeNameNoDebug("rainfallModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        rainfallModel,
        dictionary,
        (
            const fvMesh& mesh,
            const dictionary& dict
        ),
        (mesh, dict)
    );


    // Constructors

        //- Construct for given rainfallModel
        rainfallModel(const fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        rainfallModel(const rainfallModel&) = delete;


    // Selectors

        //- Select from a dictionary
        static autoPtr<rainfallModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~rainfallModel();


    // Member Functions

        //- Return the height source term
        virtual tmp<volScalarField> R0() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
