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
    Foam::fluidBlastThermo

Description
    General fluid class that contains functions to return pressure and speed
    of sound

SourceFiles
    fluidBlastThermo.C
    newfluidBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef fluidBlastThermo_H
#define fluidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidThermo.H"
#include "blastThermo.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fluidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidBlastThermo
:
    public blastThermo,
    public fluidThermo
{
protected:

        //- Pressure [Pa]
        volScalarField& p_;

        //- Dynamic viscosity [kg/m/s]
        volScalarField mu_;

        //- Speed of sound [m/s]
        volScalarField speedOfSound_;


        // protected member functions

            //- Dynamic viscosity of mixture [kg/m/s]
            volScalarField& muRef()
            {
                return mu_;
            }


public:

    //- Runtime type information
    TypeName("fluid");

    //- Construct from mesh, dictionary and phase name
    fluidBlastThermo
    (
        const fvMesh&,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    //- Initialize models
    virtual void initializeFields();


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidBlastThermo,
        dictionary,
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName,
            const word& masterName
        ),
        (mesh, dict, phaseName, masterName)
    );

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidBlastThermo,
        phase,
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName
        ),
        (mesh, dict, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<fluidBlastThermo> New
        (
            const label,
            const fvMesh&,
            const dictionary&,
            const word& phaseName = word::null
        );

    //- Destructor
    virtual ~fluidBlastThermo();

    // Member Functions

    // Access functions

        // Access to thermodynamic state variables

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p();

            //- Pressure [Pa]
            virtual const volScalarField& p() const;


        // Fields derived from thermodynamic state variables

            //- Calclate internal energy [J/kg]
            virtual tmp<volScalarField> calce(const volScalarField&) const = 0;

            //- Recompute the density
            virtual void updateRho();

            //- Recompute the density
            virtual void updateRho(const volScalarField&) = 0;

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho(const volScalarField& deltaRho)
            {}

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const
            {
                NotImplemented;
                return rho();
            }

            //- Gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const
            {
                return blastThermo::gamma();
            }

            //- Gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& T,
                const label patchi
            ) const
            {
                return blastThermo::gamma(T, patchi);
            }

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture [W/m/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField& alphat
            ) const
            {
                return blastThermo::kappaEff(alphat);
            }

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture for patch [W/m/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const
            {
                return blastThermo::kappaEff(alphat, patchi);
            }

            //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const
            {
                return blastThermo::alphaEff(alphat);
            }

            //- Effective thermal turbulent diffusivity of mixture
            //  for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const
            {
                return blastThermo::alphaEff(alphat, patchi);
            }

            using blastThermo::gamma;
            using blastThermo::he;

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Calculate pressure for cell [Pa]
            virtual scalar cellpRhoT(const label celli) const = 0;

            //- Return the speed of sound
            virtual const volScalarField& speedOfSound() const;

            //- Access the speed of sound
            virtual volScalarField& speedOfSound();

            //- Derivative of pressure w.r.t density for cel
            virtual scalar celldpdRho(const label celli) const = 0;

            //- Derivative of pressure w.r.t internal energy for cell
            virtual scalar celldpde(const label celli) const = 0;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual tmp<volScalarField> mu() const;

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> mu(const label patchi) const;

            using fluidThermo::nu;

            //- Kinematic viscosity of [m^2/s]
            virtual tmp<volScalarField> nu() const
            {
                return fluidThermo::nu();
            }

            //- Kinematic viscosity of mixture for cell [m^2/s]
            virtual scalar cellnu(const label celli) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
