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
    Foam::multicomponentBlastThermo

Description
    Class to handle species transport

SourceFiles
    multicomponentBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentBlastThermo_H
#define multicomponentBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicBlastThermo.H"
#include "PtrListDictionary.H"
#include "speciesTable.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multicomponentBlastThermo Declaration
\*---------------------------------------------------------------------------*/

class multicomponentBlastThermo
:
    private integrationSystem
{
//- Private data

    // Const reference to the mesh
    const fvMesh& mesh_;

    //- Names of phase
    word phaseName_;

    //- Name of master thermo class
    word masterName_;

protected:
//- Protected data

    //- Name of species
    speciesTable species_;

    //- Mass fractions
    PtrList<volScalarField> Ys_;

    //- Mass transfer rates passed from an external source
    PtrListDictionary<volScalarField> massTransferRates_;

    //- Index of inert specie
    label inertIndex_;

    //- List of active species
    List<bool> active_;

    //- Name of flux
    UautoPtr<const surfaceScalarField> alphaRhoPhiPtr_;

    //- Name of mass
    UautoPtr<const volScalarField> alphaRhoPtr_;


public:

    // Constructor
    multicomponentBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    // Constructor
    multicomponentBlastThermo
    (
        const speciesTable& species,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );


    //- Destructor
    virtual ~multicomponentBlastThermo();


    //- ODE functions

        //- Initialize models
        virtual void initializeModels();

        //- Solve sub-step stepi
        virtual void solve();

        //- Add an external change to a thermodynamic field if present
        //  i.e. mass transfer
        virtual void addDelta(const word& name, const volScalarField& delta);


    // Access

        //- Return the table of species
        virtual speciesTable species() const
        {
            return species_;
        }

        //- Does the mixture include this specie?
        virtual bool contains(const word& specieName) const
        {
            return species_[specieName];
        }

        //- Return true for active species
        inline bool active(label speciei) const;

        //- Return the bool list of active species
        inline const List<bool>& active() const;

        //- Set speciei active
        inline void setActive(label speciei);

        //- Set speciei inactive
        inline void setInactive(label speciei);

        //- Return the mass-fraction fields
        inline PtrList<volScalarField>& Ys();

        //- Return the const mass-fraction fields
        inline const PtrList<volScalarField>& Ys() const;

        //- Return the mass-fraction field for a specie given by index
        inline volScalarField& Y(const label i);

        //- Return the const mass-fraction field for a specie given by index
        inline const volScalarField& Y(const label i) const;

        //- Return the mass-fraction field for a specie given by name
        inline volScalarField& Y(const word& specieName);

        //- Return the const mass-fraction field for a specie given by name
        inline const volScalarField& Y(const word& specieName) const;

        //- Return the name of the thermo physics
        virtual word thermoName() const = 0;


        //- Return species properties

            //- Return molecular weight
            virtual tmp<volScalarField> specieW(const word& s) const = 0;

            //- Return molecular weight
            virtual scalar specieWi(const word& s, const label celli) const = 0;

            //- Mie Gruniesen coefficient
            virtual tmp<volScalarField> specieGamma(const word& s) const = 0;

            //- Mie Gruniesen coefficient for a cell
            virtual scalar specieGammai(const word& s, const label celli) const = 0;

            //- Heat capacity at constant pressure of mixture [J/kg/K]
            virtual tmp<volScalarField> specieCp(const word& s) const = 0;

            //- Heat capacity at constant pressure of mixture for a cell [J/kg/K]
            virtual scalar specieCpi(const word& s, const label celli) const = 0;

            //- Heat capacity at constant volume of mixture [J/kg/K]
            virtual tmp<volScalarField> specieCv(const word& s) const = 0;

            //- Heat capacity at constant volume of mixture for a cell [J/kg/K]
            virtual scalar specieCvi(const word& s, const label celli) const = 0;

            //- Specific heat ratio
            virtual tmp<volScalarField> specieCpByCv(const word& s) const = 0;

            //- Heat of formation
            virtual tmp<volScalarField> specieHf(const word& s) const = 0;

            //- Heat of formation
            virtual tmp<volScalarField> specieFlameT(const word& s) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentBlastThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
