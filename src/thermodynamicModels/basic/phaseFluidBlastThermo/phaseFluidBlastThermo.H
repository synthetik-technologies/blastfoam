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
    Foam::phaseFluidBlastThermo

Description
    General fluid class that contains functions to return pressure and speed
    of sound

SourceFiles
    phaseFluidBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef phaseFluidBlastThermo_H
#define phaseFluidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicBlastThermo.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseFluidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

class phaseFluidBlastThermo
:
    public basicBlastThermo
{
protected:

    // Protected data

        //- Pressure [Pa]
        const volScalarField& p_;


public:

    //- Runtime type information
    TypeName("phaseFluid");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseFluidBlastThermo,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const word& phaseName,
                const word& masterName
            ),
            (mesh, dict, phaseName, masterName)
        );


    // Constructor

        //- Construct with input field
        phaseFluidBlastThermo
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName,
            const word& masterName = word::null
        );

    //- Initialize models
    virtual void initializeModels()
    {}

    // Selectors

        static autoPtr<phaseFluidBlastThermo> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName,
            const word& masterName = word::null
        );

    //- Destructor
    virtual ~phaseFluidBlastThermo();


    //- Thermodynamic and transport functions

        // Access to derived thermodynamic variables

            //- Calculate density
            virtual void correct();

            //- Calculate density
            virtual void updateRho(const volScalarField&) = 0;

            //- Calculate internal energy given a pressure
            virtual tmp<volScalarField> calce(const volScalarField&) const = 0;

            //- Pressure [Pa]
            const volScalarField& p() const;

            //- Return the speed of sound
            virtual tmp<volScalarField> pRhoT() const = 0;

            //- Calculate pressure for cell [Pa]
            virtual scalar cellpRhoT(const label celli) const = 0;

            //- Derivative of pressure w.r.t density for cell
            virtual scalar celldpdRho(const label celli) const = 0;

            //- Derivative of pressure w.r.t internal energy for cell
            virtual scalar celldpde(const label celli) const = 0;

            //- Return the speed of sound
            virtual tmp<volScalarField> speedOfSound() const = 0;

            //- Return the speed of sound for patchi
            virtual tmp<scalarField> speedOfSound(const label patchi) const = 0;


        // Fields derived from thermodynamic state variables

            //- Gamma = dpde []
            virtual tmp<volScalarField> Gamma() const = 0;

            //- Gamma = dpde for patch []
            virtual tmp<scalarField> Gamma(const label patchi) const = 0;

            //- Gamma = dpde for cell []
            virtual scalar cellGamma(const label celli) const = 0;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual tmp<volScalarField> calcMu() const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
