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
    Foam::blastRadiationModel

Description
    Top level model for radiation modelling

SourceFiles
    blastRadiationModel.C
    blastRadiationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef blastRadiationModel_H
#define blastRadiationModel_H

#include "radiationModel.H"
#include "runTimeSelectionTables.H"
#include "radiationODE.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace radiationModels
{
class blastAbsorptionEmissionModel;
}


/*---------------------------------------------------------------------------*\
                       Class blastRadiationModel Declaration
\*---------------------------------------------------------------------------*/

class blastRadiationModel
:
    public radiationModel
{
protected:

    // Protected data

        //- Blast absorption/emission model (type cast)
            UautoPtr<const radiationModels::blastAbsorptionEmissionModel>
                bAbsorptionEmission_;

        //- Radiation ode
        radiationODE radODE_;

private:

    // Private Member Functions

        //- Initialise
        void initialise();


public:

    //- Runtime type information
    TypeName("blastRadiationModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            blastRadiationModel,
            T,
            (
                const volScalarField& T
            ),
            (T)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            blastRadiationModel,
            dictionary,
            (
                const dictionary& dict,
                const volScalarField& T
            ),
            (dict, T)
        );


    // Constructors

        //- Null constructor
        blastRadiationModel(const volScalarField& T);

        //- Construct from components
        blastRadiationModel(const word& type, const volScalarField& T);

        //- Construct from components
        blastRadiationModel
        (
            const word& type,
            const dictionary& dict,
            const volScalarField& T
        );

        //- Disallow default bitwise copy construction
        blastRadiationModel(const blastRadiationModel&) = delete;


    // Selectors

        //- Return a reference to the selected radiation model
        static autoPtr<blastRadiationModel> New(const volScalarField& T);

        //- Return a reference to the selected radiation model
        static autoPtr<blastRadiationModel> New
        (
            const dictionary& dict,
            const volScalarField& T
        );


    //- Destructor
    virtual ~blastRadiationModel();


    // Member Functions

        // Access

            //- Source term component (for power of T^4)
            virtual scalar Rpi(const label celli) const = 0;

            //- Source term component (constant)
            virtual scalar Rui(const label celli) const = 0;

            //- Source term component (constant)
            virtual tmp<volScalarField> calcRhoE
            (
                const dimensionedScalar& dt,
                const volScalarField& rhoE,
                const volScalarField& rho,
                const volScalarField& e,
                const volScalarField& Cv
            );

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const blastRadiationModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addToBlastRadiationRunTimeSelectionTables(model)                       \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        blastRadiationModel,                                                   \
        model,                                                                 \
        dictionary                                                             \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        blastRadiationModel,                                                   \
        model,                                                                 \
        T                                                                      \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
