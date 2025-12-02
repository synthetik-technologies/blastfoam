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
#include "physicalProperties.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fluidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidBlastThermo
:
    public physicalProperties,
    public blastThermo,
    virtual public fluidThermo
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

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            volScalarField& pRef()
            {
                return p_;
            }

            //- Access the speed of sound
            volScalarField& speedOfSoundRef()
            {
                return speedOfSound_;
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
        const word& masterName = word::null,
        const bool requireRho = true
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
            const fvMesh&,
            const dictionary&,
            const word& thermoType,
            const word& phaseName = word::null
        );

        //- Standard selection based on fvMesh
        static autoPtr<fluidBlastThermo> New
        (
            const fvMesh&,
            const word& thermoType,
            const word& phaseName = word::null
        );

    //- Destructor
    virtual ~fluidBlastThermo();

    // Member Functions

    // Access functions

        //- Properties dictionary
        virtual const IOdictionary& properties() const
        {
            return *this;
        }

        //- Non-const access the properties dictionary
        virtual IOdictionary& properties()
        {
            return *this;
        }
        virtual void clear()
        {
            return blastThermo::clear();
        }

        // Access to thermodynamic state variables

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p();

            //- Pressure [Pa]
            virtual const volScalarField& p() const;

            virtual tmp<volScalarField> rho() const
            {
                return rho_;
            }

            //- Rename and return the thermodynamic density field [kg/m^3]
            //  This is used by solvers which create a separate continuity rho
            virtual tmp<volScalarField> renameRho();


        // Fields derived from thermodynamic state variables

            //- Calculate internal energy [J/kg]
            virtual tmp<volScalarField> calce(const volScalarField&) const = 0;

            //- Calculate internal energy for a cell
            virtual scalar calcCelle
            (
                const scalar p,
                const label celli
            ) const = 0;

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
            virtual scalar cellpRhoT
            (
                const label celli,
                const bool limit = true
            ) const = 0;

            //- Calculate pressure for patch and face [Pa]
            virtual scalar patchFacepRhoT
            (
                const label patchi,
                const label facei,
                const bool limit = true
            ) const = 0;

            //- Calculate the pressure
            virtual tmp<volScalarField> pRhoT() const;

            //- Return the speed of sound
            virtual const volScalarField& speedOfSound() const;

            //- Access the speed of sound
            virtual volScalarField& speedOfSound();

            //- Calculate the speed of sound in a cell
            virtual scalar cellSpeedOfSound
            (
                const scalar p,
                const label celli
            ) const = 0;

            //- Derivative of pressure w.r.t density for cel
            virtual scalar celldpdRho(const label celli) const = 0;

            //- Derivative of pressure w.r.t internal energy for cell
            virtual scalar celldpde(const label celli) const = 0;

            //- Derivative of pressure w.r.t temperature for cell
            virtual scalar celldpdT(const label celli) const = 0;

            //- Is the equation of state temperature based
            virtual bool temperatureBased() const = 0;


        // Access to transport state variables

            //- Return if inviscid
            virtual bool inviscid() const = 0;

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual const volScalarField& mu() const;

            //- Kinematic viscosity of [m^2/s]
            using fluidThermo::nu;

            //- Kinematic viscosity of [m^2/s]
            virtual tmp<volScalarField> nu() const
            {
                return fluidThermo::nu();
            }

            //- Kinematic viscosity of mixture for cell [m^2/s]
            virtual scalar cellnu(const label celli) const;


        //- Read thermophysical properties dictionary
        virtual bool read()
        {
            return blastThermo::read();
        }
        virtual bool read(const dictionary& dict)
        {
            return blastThermo::read(dict);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
