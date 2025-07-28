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
    Foam::solidBlastThermo

Description
    General solid class

SourceFiles
    solidBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidBlastThermo_H
#define solidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidThermo.H"
#include "blastThermo.H"
#include "physicalProperties.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class solidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

class solidBlastThermo
:
    public physicalProperties,
    public blastThermo,
    public solidThermo
{
protected:

    //- Anisotropic thermal conductivity [W/m/K]
    autoPtr<volVectorField> KappaPtr_;

    //- Thermal conductivity for enthalpy of mixture [kg/m/s]
    volVectorField& KappaRef()
    {
        return KappaPtr_();
    }


public:

    //- Runtime type information
    TypeName("solid");

    //- Construct from mesh, dictionary and phase name
    solidBlastThermo
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
        solidBlastThermo,
        dictionary,
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName,
            const word& masterName
        ),
        (mesh, dict, phaseName, masterName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<solidBlastThermo> New
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName = word::null
        );

        //- Standard selection based on fvMesh
        static autoPtr<solidBlastThermo> New
        (
            const fvMesh&,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~solidBlastThermo();

    // Member Functions

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

        //- The dpdt term should not be included in the enthalpy equation
        virtual Switch dpdt() const
        {
            return false;
        }

        // Fields derived from thermodynamic state variables

            //- Calculate internal energy given a pressure
            virtual tmp<volScalarField> calce() const = 0;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho()
            {
                return blastThermo::rho();
            }

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const
            {
                return blastThermo::rho();
            }

            //- Density for a patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const
            {
                return blastThermo::rho(patchi);
            }

            //- Old-time density [kg/m^3]
            virtual tmp<volScalarField> rho0() const
            {
                return blastThermo::rho0();
            }

            using blastThermo::rho;

            //- Kinematic viscosity [m^2/s]
            virtual tmp<volScalarField> nu() const;


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


        // Access to transport state variables

            //- Anisotropic thermal conductivity for patch
            //  in the local coordinate system [W/m/K]
            virtual const volVectorField& Kappa() const
            {
                return KappaPtr_();
            }


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
