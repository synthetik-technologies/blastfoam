/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technologies
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
    Foam::atmosphereModel

Description
    Base class for atmosphere model

SourceFiles
    atmosphereModel.C
    newAtmosphereModel.C

\*---------------------------------------------------------------------------*/

#ifndef atmosphereModel_H
#define atmosphereModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "uniformDimensionedFields.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class atmosphereModel Declaration
\*---------------------------------------------------------------------------*/

class atmosphereModel
{
protected:

    //- Reference to atmosphereProperties
    const dictionary& dict_;

    //- Const reference to mesh
    const fvMesh& mesh_;

    //- Gravitational acceleration
    uniformDimensionedVectorField g_;

    //- Normal vector for gravitational acceleration
    vector normal_;

    //- Elevation of lowest point in mesh
    dimensionedScalar groundElevation_;

    //- Height field
    volScalarField h_;


public:

    //- Runtime type information
    TypeName("atmosphereModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            atmosphereModel,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );

    // Constructor
    atmosphereModel(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~atmosphereModel();


    // Selectors

        static autoPtr<atmosphereModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        virtual void createAtmosphere
        (
            volScalarField& p,
            volScalarField& rho
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
