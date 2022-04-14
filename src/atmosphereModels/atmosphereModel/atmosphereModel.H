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
    atmosphereModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef atmosphereModel_H
#define atmosphereModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "uniformDimensionedFields.H"
#include "fluidBlastThermo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class atmosphereModel Declaration
\*---------------------------------------------------------------------------*/

class atmosphereModel
{
protected:

    //- Reference to atmosphereProperties
    dictionary dict_;

    //- Const reference to mesh
    const fvMesh& mesh_;

    //- Gravitational acceleration
    uniformDimensionedVectorField g_;

    //- Normal vector for gravitational acceleration
    vector normal_;

    //- Elevation of lowest point in mesh
    uniformDimensionedScalarField baseElevation_;

    //- Height field
    volScalarField h_;

    //- Saved zoneID to set
    const label zoneID_;

    //- Stabilize fields
    void hydrostaticInitialisation
    (
        fluidBlastThermo& thermo
    ) const;

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
                const dictionary& dict,
                const label zoneID
            ),
            (mesh, dict, zoneID)
        );

    // Constructor
    atmosphereModel(const fvMesh&, const dictionary&, const label);


    //- Destructor
    virtual ~atmosphereModel();


    // Selectors

        // Construct from dictionary
        static autoPtr<atmosphereModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const label zoneID = -1
        );

        //- COnstruct given a model type and dictionary
        static autoPtr<atmosphereModel> New
        (
            const word& atmosphereModelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const label zoneID = -1
        );


    // Member Functions

        virtual void createAtmosphere(fluidBlastThermo& thermo) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
