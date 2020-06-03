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
    Foam::fluidThermoModel

Description
    General class for a fluid/solid

SourceFiles
    fluidThermoModel.C

\*---------------------------------------------------------------------------*/

#ifndef fluidThermoModel_H
#define fluidThermoModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicThermoModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "wordIOList.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluidThermoModel Declaration
\*---------------------------------------------------------------------------*/

class fluidThermoModel
:
    public basicThermoModel
{
protected:

    // Protected data

        //- Saved viscosity field
        volScalarField mu_;

        //- Residual volume fraction
        dimensionedScalar residualAlpha_;

        //- Residual density
        dimensionedScalar residualRho_;

        //- Is there viscosity
        bool viscous_;

        //- Initialize internal energy
        void initialize();

public:

    //- Runtime type information
    TypeName("fluidThermo");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidThermoModel,
            basic,
            (
                const word& name,
                volScalarField& p,
                volScalarField& rho,
                volScalarField& e,
                volScalarField& T,
                const dictionary& dict,
                const bool master
            ),
            (name, p, rho, e, T, dict, master)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidThermoModel,
            detonating,
            (
                const word& name,
                volScalarField& p,
                volScalarField& rho,
                volScalarField& e,
                volScalarField& T,
                const dictionary& dict,
                const bool master
            ),
            (name, p, rho, e, T, dict, master)
        );

    // Constructor

        //- Construct with input field
        fluidThermoModel
        (
            const word& phaseName,
            volScalarField& p,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const bool master = true
        );

        //- Construct and add fields to db
        fluidThermoModel
        (
            const word& phaseName,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool master = true
        );


    // Selectors

        static autoPtr<fluidThermoModel> NewBasic
        (
            const word& name,
            volScalarField& p,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const bool master = true
        );

        static autoPtr<fluidThermoModel> NewDetonating
        (
            const word& name,
            volScalarField& p,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const bool master = true
        );

        static autoPtr<fluidThermoModel> New
        (
            const word& name,
            volScalarField& p,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const bool master = true
        );

    //- Destructor
    virtual ~fluidThermoModel();

    // Access functions

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

        //- Return the speed of sound
        virtual tmp<volScalarField> speedOfSound() const = 0;

        //- Return the speed of sound for patchi
        virtual tmp<scalarField> speedOfSound(const label patchi) const = 0;

        //- Non-const access to thermo model phasei
        virtual fluidThermoModel& thermo(const label)
        {
            return *this;
        }

        //- Const access to thermo model phasei
        virtual const fluidThermoModel& thermo(const label) const
        {
            return *this;
        }

        const dimensionedScalar& residualAlpha() const
        {
            return residualAlpha_;
        }

        const dimensionedScalar& residualRho() const
        {
            return residualRho_;
        }


    //- Thermodynamic and transport functions

        //- Calculate thermodynamic pressure
        virtual tmp<volScalarField> calcP() const = 0;

        //- Calculate internal energy from pressure and density
        virtual tmp<volScalarField> calce() const = 0;

        //- Dynamic viscosity of mixture [kg/m/s]
        virtual tmp<volScalarField> mu() const;

        //- Dynamic viscosity of mixture for patch [kg/m/s]
        virtual tmp<scalarField> mu(const label patchi) const;

        //- Kinematic viscosity of mixture [m^2/s]
        virtual tmp<volScalarField> nu() const;

        //- Kinematic viscosity of mixture for patch [m^2/s]
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Kinematic viscosity of mixture for cell [m^2/s]
        virtual scalar nui(const label celli) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
