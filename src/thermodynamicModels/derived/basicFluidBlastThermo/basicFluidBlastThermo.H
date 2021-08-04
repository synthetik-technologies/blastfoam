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
    Foam::basicFluidBlastThermo

Description
    Non reacting, single phase fluid thermo model

SourceFiles
    basicFluidBlastThermo.C
    basicFluidBlastThermos.C

\*---------------------------------------------------------------------------*/

#ifndef basicFluidBlastThermo_H
#define basicFluidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidBlastThermo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class basicFluidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class basicFluidBlastThermo
:
    public Thermo
{
public:

    //- Runtime type information
    TypeName
    (
        word("basic<fluid<" + Thermo::thermoType::typeName() + ">>").c_str()
    );

    // Constructor
    basicFluidBlastThermo
    (
        const word& phaseName,
        volScalarField& rho,
        volScalarField& e,
        volScalarField& T,
        const dictionary& dict,
        const word& masterName = word::null
    );


    //- Destructor
    virtual ~basicFluidBlastThermo();


    virtual word thermoName() const ;
    // Member Functions

        //- Pre update
        virtual void update()
        {}

        //- Solve sub-step stepi
        virtual void solve()
        {}

        //- Remove stored fields
        virtual void clearODEFields()
        {}

        //- Post update
        virtual void postUpdate()
        {}

        //- Correct fields
        virtual void correct();

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;

        //- Return the speed of sound
        virtual tmp<volScalarField> speedOfSound() const;

        //- Return the speed of sound for patchi
        virtual tmp<scalarField> speedOfSound(const label patchi) const;

        //- Gamma = dpde []
        virtual tmp<volScalarField> Gamma() const;

        //- Gamma = dpde for patch []
        virtual tmp<scalarField> Gamma(const label patchi) const;

        //- Gamma = dpde for cell []
        virtual scalar Gammai(const label celli) const;


    //- Thermodynamic and transport functions

        //- Dynamic viscosity of mixture [kg/m/s]
        virtual tmp<volScalarField> mu() const;

        //- Dynamic viscosity of mixture for patch [kg/m/s]
        virtual tmp<scalarField> mu(const label patchi) const;

        //- Calculate density
        virtual void updateRho(const volScalarField& p);

        //- Calculate thermodynamic pressure
        virtual tmp<volScalarField> pRhoT() const;

        //- Calculate internal energy from pressure and density
        virtual tmp<volScalarField> calce(const volScalarField& p) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicFluidBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
