/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::blastCombustionModel

Description
    Base class for combustion models

SourceFiles
    blastCombustionModel.C

\*---------------------------------------------------------------------------*/

#ifndef blastCombustionModel_H
#define blastCombustionModel_H

#include "IOdictionary.H"
#include "compressibleMomentumTransportModels.H"
#include "multicomponentBlastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class blastCombustionModel Declaration
\*---------------------------------------------------------------------------*/

class blastCombustionModel
:
    public IOdictionary
{
    // Private Member Functions

        //- Construct the base IO object
        IOobject createIOobject
        (
            const multicomponentBlastThermo& thermo,
            const word& combustionProperties
        ) const;


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the thermo
        const multicomponentBlastThermo& thermo_;

        //- Dictionary of the model
        dictionary coeffs_;

        //- Model type
        const word modelType_;


public:

    //- Runtime type information
    TypeName("combustionModel");

    //- Default combustionProperties dictionary name
    static const word combustionPropertiesName;


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        blastCombustionModel,
        fluid,
        (
            const word& modelType,
            const multicomponentBlastThermo& thermo,
            const bool isFluid,
            const word& combustionProperties
        ),
        (modelType, thermo, isFluid, combustionProperties)
    );
    declareRunTimeSelectionTable
    (
        autoPtr,
        blastCombustionModel,
        solid,
        (
            const word& modelType,
            const multicomponentBlastThermo& thermo,
            const bool isFluid,
            const word& combustionProperties
        ),
        (modelType, thermo, isFluid, combustionProperties)
    );


    // Constructors

        //- Construct from components
        blastCombustionModel
        (
            const word& modelType,
            const multicomponentBlastThermo& thermo,
            const word& combustionProperties = combustionPropertiesName
        );

        //- Disallow default bitwise copy construction
        blastCombustionModel(const blastCombustionModel&) = delete;


    // Selectors

        //- Select using thermo and turbulence
        static autoPtr<blastCombustionModel> NewFluid
        (
            const multicomponentBlastThermo& thermo,
            const word& combustionProperties = combustionPropertiesName,
            const bool required = true
        );
        //- Select using thermo and turbulence
        static autoPtr<blastCombustionModel> NewSolid
        (
            const multicomponentBlastThermo& thermo,
            const word& combustionProperties = combustionPropertiesName,
            const bool required = true
        );


    //- Destructor
    virtual ~blastCombustionModel();


    // Member Functions

        //- Return const access to the mesh database
        inline const fvMesh& mesh() const;

        //- Return const access to the thermo
        inline const multicomponentBlastThermo& thermo() const;

        //- Return access to turbulence
        inline const compressibleMomentumTransportModel& turbulence() const;

        //- Return const access to rho
        inline const volScalarField& rho() const;

        //- Return const access to phi
        inline tmp<surfaceScalarField> phi() const;

        //- Return const dictionary of the model
        inline const dictionary& coeffs() const;

        //- Correct combustion rate
        virtual void correct() = 0;

        //- Consumption rate, i.e. source term for specie equation
        virtual tmp<volScalarField::Internal> R(const label speciei) const = 0;

        //- Consumption rate matrix, i.e. source term for specie equation
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const = 0;

        //- Heat release rate [kg/m/s^3]
        virtual tmp<volScalarField> Qdot() const = 0;

        //- Update properties from given dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const blastCombustionModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "blastCombustionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
