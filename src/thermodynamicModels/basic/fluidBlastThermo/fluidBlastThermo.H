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

#include "basicBlastThermo.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "wordIOList.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidBlastThermo
:
    public basicBlastThermo
{
protected:

    // Protected data

        //- Saved viscosity field
        volScalarField mu_;

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
            fluidBlastThermo,
            basic,
            (
                const word& name,
                volScalarField& p,
                volScalarField& rho,
                volScalarField& e,
                volScalarField& T,
                const dictionary& dict,
                const bool master,
                const word& masterName
            ),
            (name, p, rho, e, T, dict, master, masterName)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidBlastThermo,
            detonating,
            (
                const word& name,
                volScalarField& p,
                volScalarField& rho,
                volScalarField& e,
                volScalarField& T,
                const dictionary& dict,
                const bool master,
                const word& masterName
            ),
            (name, p, rho, e, T, dict, master, masterName)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidBlastThermo,
            multicomponent,
            (
                const word& name,
                volScalarField& p,
                volScalarField& rho,
                volScalarField& e,
                volScalarField& T,
                const dictionary& dict,
                const bool master,
                const word& masterName
            ),
            (name, p, rho, e, T, dict, master, masterName)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidBlastThermo,
            reacting,
            (
                const word& name,
                volScalarField& p,
                volScalarField& rho,
                volScalarField& e,
                volScalarField& T,
                const dictionary& dict,
                const bool master,
                const word& masterName
            ),
            (name, p, rho, e, T, dict, master, masterName)
        );

    // Constructor

        //- Construct with input field
        fluidBlastThermo
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
        fluidBlastThermo
        (
            const word& phaseName,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );


    // Selectors

        static autoPtr<fluidBlastThermo> NewBasic
        (
            const word& name,
            volScalarField& p,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );

        static autoPtr<fluidBlastThermo> NewDetonating
        (
            const word& name,
            volScalarField& p,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );

        static autoPtr<fluidBlastThermo> NewMulticomponent
        (
            const word& name,
            volScalarField& p,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );

        static autoPtr<fluidBlastThermo> NewReacting
        (
            const word& name,
            volScalarField& p,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );

        static autoPtr<fluidBlastThermo> New
        (
            const word& name,
            volScalarField& p,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const bool master = true,
            const word& masterName = word::null
        );

    //- Destructor
    virtual ~fluidBlastThermo();

    // Access functions

        //- Return energy source
        virtual tmp<volScalarField> ESource() const = 0;

        //- Return the speed of sound
        virtual tmp<volScalarField> speedOfSound() const = 0;

        //- Return the speed of sound for patchi
        virtual tmp<scalarField> speedOfSound(const label patchi) const = 0;

        //- Non-const access to thermo model phasei
        virtual fluidBlastThermo& thermo(const label)
        {
            return *this;
        }

        //- Const access to thermo model phasei
        virtual const fluidBlastThermo& thermo(const label) const
        {
            return *this;
        }


    //- Thermodynamic and transport functions

        //- Calculate thermodynamic pressure
        virtual tmp<volScalarField> calcP() const;

        //- Calculate thermodynamic pressure for a patch
        virtual tmp<scalarField> calcP(const label patchi) const = 0;

        //- Calculate thermodynamic pressure
        virtual scalar calcPi(const label celli) const = 0;

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
