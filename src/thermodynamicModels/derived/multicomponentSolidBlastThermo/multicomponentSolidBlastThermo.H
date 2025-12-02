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
    Foam::multicomponentSolidBlastThermo

Description
    Multicomponent fluid model with species transport

SourceFiles
    multicomponentSolidBlastThermo.C
    multicomponentSolidBlastThermos.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentSolidBlastThermo_H
#define multicomponentSolidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentBlastThermo.H"
#include "mixtureBlastThermo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class multicomponentSolidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class multicomponentSolidBlastThermo
:
    public Thermo
{
    //- Add contribution to mixture temperature
    virtual void calculate();


public:

    //- Runtime type information
    TypeName("multicomponent");

    // Constructor
    multicomponentSolidBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    // Constructor
    multicomponentSolidBlastThermo
    (
        const HashPtrTable<Thermo, word, string::hash>& thermoData,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    //- Destructor
    virtual ~multicomponentSolidBlastThermo();


    // Member Functions

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;

        //- Calculate internal energy
        virtual tmp<volScalarField> calce() const;

        //- Correct thermodynamic fields
        virtual void correct();

        //- Pre update step
        virtual void update()
        {
            multicomponentBlastThermo::update();
        }

        //- Solve sub-step stepi
        virtual void solve()
        {
            multicomponentBlastThermo::solve();
        }

        //- Explicit post update
        virtual void postExplicit()
        {
            multicomponentBlastThermo::postExplicit();
        }

        //- Implicit post update
        virtual void postImplicit()
        {
            multicomponentBlastThermo::postImplicit();
        }

        //- Store explicit deltas
        virtual void storeExplicit()
        {
            multicomponentBlastThermo::storeExplicit();
        }

        //- Clear stored data
        virtual void clear()
        {
            multicomponentBlastThermo::clear();
        }

    // Specie thermo properties

            //- Return the pressure for cell [Pa]
            virtual scalar cellp(const label celli) const;

            //- Pressure [kg/m/s^2]
            virtual scalar pi
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const;

            //- Density [kg/m^3]
            virtual tmp<volScalarField> pi
            (
                const label speciei,
                const volScalarField& rho,
                const volScalarField& e,
                const volScalarField& T
            ) const;

            //- Return derivative of pressure w.r.t. density
            virtual scalar dpdRhoi
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const;

            //- Return derivative of pressure w.r.t. Temperature
            virtual scalar dpdTi
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const;


        // Mass Transfer

            //- Does the mixture include this specie?
            inline bool containsSpecie(const word& specieName) const
            {
                return multicomponentBlastThermo::containsSpecie(specieName);
            }
            //- Add an external change to a thermodynamic field if present
            //  i.e. mass transfer
            virtual void addDelta(const word& name, tmp<volScalarField>& delta)
            {
                multicomponentBlastThermo::addDelta(name, delta);
            }

            //- Add an external change to a thermodynamic field if present
            //  i.e. mass transfer
            virtual void addDelta(const word& name, const volScalarField::Internal& delta)
            {
                multicomponentBlastThermo::addDelta(name, delta);
            }

            //- Add an external change to a thermodynamic field if present
            //  i.e. mass transfer
            virtual void addSource(const word& name, tmp<fvScalarMatrix>& source)
            {
                multicomponentBlastThermo::addSource(name, source);
            }

    //- Thermodynamic and transport functions

        //- Calculate density
        virtual void updateRho();

        //- Dynamic viscosity [kg/m/s]
        virtual scalar mui
        (
            const label speciei,
            const scalar p,
            const scalar T
        ) const;

        //- Dynamic viscosity [kg/m/s]
        virtual scalar mui
        (
            const label speciei,
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Dynamic viscosity [kg/m/s]
        virtual tmp<volScalarField> mui
        (
            const label speciei,
            const volScalarField& p,
            const volScalarField& T
        ) const;

        //- Dynamic viscosity [kg/m/s]
        virtual tmp<volScalarField> mui
        (
            const label speciei,
            const volScalarField& rho,
            const volScalarField& e,
            const volScalarField& T
        ) const;

        //- Return true if thermal conductivity is isotropic
        virtual bool isotropic() const
        {
            return Thermo::thermoType::isotropic;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multicomponentSolidBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
