/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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
    Foam::solidThermoModel

Description
    Template class to define a runtime selectable fluidThermo model

SourceFiles
    solidThermoModel.C

\*---------------------------------------------------------------------------*/

#ifndef solidThermoModel_H
#define solidThermoModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicThermoModel.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class basicFluidThermo Declaration
\*---------------------------------------------------------------------------*/

class solidThermoModel
:
    public basicThermoModel
{

public:

    //- Runtime type information
    TypeName("solidThermoModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        solidThermoModel,
        basicSolid,
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool master,
            const word& masterName
        ),
        (name, mesh, dict, master, masterName)
    );

    declareRunTimeSelectionTable
    (
        autoPtr,
        solidThermoModel,
        detonatingSolid,
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool master,
            const word& masterName
        ),
        (name, mesh, dict, master, masterName)
    );

    // Constructor

        //- Construct with input field
        solidThermoModel
        (
            const word& phaseName,
            volScalarField& p,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );

        //- Construct and add fields to db
        solidThermoModel
        (
            const word& phaseName,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );

    //- Selector
        static autoPtr<solidThermoModel> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );

        static autoPtr<solidThermoModel> NewBasic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );

        static autoPtr<solidThermoModel> NewDetonating
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );

    //- Destructor
    virtual ~solidThermoModel();


        //- Solve sub-step stepi
        virtual void solve
        (
            const label stepi,
            const scalarList& ai,
            const scalarList& bi
        ) = 0;

        //- Set old lists and fluxes (initialization of fields)
        virtual void setODEFields
        (
            const label nSteps,
            const labelList& oldIs,
            const label& nOld,
            const labelList& deltaIs,
            const label nDelta
        ) = 0;

        //- Remove stored fields
        virtual void clearODEFields() = 0;

        //- Return energy source
        virtual tmp<volScalarField> ESource() const = 0;

        //- Thermal conductivity [W/m/K]
        virtual tmp<volVectorField> Kappa() const = 0;

        //- Anisotropic thermal conductivity [W/m/K]
        virtual tmp<vectorField> Kappa
        (
            const label patchi
        ) const = 0;

        //- Return true if thermal conductivity is isotropic
        virtual bool isotropic() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
