/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2022
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
    const dictionary& dict_;

    //- Const reference to mesh
    const fvMesh& mesh_;

    //- Gravitational acceleration
    uniformDimensionedVectorField g_;

    //- Elevation of lowest point in mesh
    uniformDimensionedScalarField hRef_;

    //- Saved zoneID to set
    const word zoneName_;

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
                const word& zoneName
            ),
            (mesh, dict, zoneName)
        );

    // Constructor
    atmosphereModel(const fvMesh&, const dictionary&, const word&);


    //- Destructor
    virtual ~atmosphereModel();


    // Selectors

        // Construct from dictionary
        static autoPtr<atmosphereModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& zoneName = word::null
        );

        //- COnstruct given a model type and dictionary
        static autoPtr<atmosphereModel> New
        (
            const word& atmosphereModelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& zoneName = word::null
        );


    // Member Functions

        //- Return the coefficient dictionary
        const dictionary& dict() const
        {
            return dict_;
        }

        virtual void createAtmosphere(fluidBlastThermo& thermo) const = 0;

        virtual bool needReference() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
