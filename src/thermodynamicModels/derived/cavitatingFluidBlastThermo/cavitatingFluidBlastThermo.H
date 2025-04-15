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
    Foam::cavitatingFluidBlastThermo

Description
    Detonating fluid model with options to include activation and afterburn
    models

SourceFiles
    cavitatingFluidBlastThermo.C
    cavitatingFluidBlastThermos.C

\*---------------------------------------------------------------------------*/

#ifndef cavitatingFluidBlastThermo_H
#define cavitatingFluidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "blendedBlastThermo.H"
#include "saturationPressureModel.H"
#include "lookupTable1D.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cavitatingFluidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class cavitatingFluidBlastThermo
:
    public Thermo
{
// Private data

    //- Vapor fraction
    volScalarField x_;

    //- Saturation pressure [Pa] as a function of temperature
    autoPtr<saturationPressureModel> pSat_;

    //- Heat of vaporization
    dimensionedScalar hv_;

    //- Saturation density lookup table for vapor phase
    lookupTable1D<scalar> rhoSatl_;

    //- Saturation density lookup table for liquid phase
    lookupTable1D<scalar> rhoSatv_;

    //- ODE variables
    tmp<volScalarField> alphaRhoXOld_;
    tmp<volScalarField> deltaAlphaRhoX_;
    tmp<volScalarField> dxdt_;


    //- Correct thermodynamic fields
    void calculate();

    //- Return the blending field for patchi
    virtual tmp<scalarField> x(const label patchi) const
    {
        return x_.boundaryField()[patchi];
    }

    //- Return the blending field for celli
    virtual scalar cellx(const label celli) const
    {
        return x_[celli];
    }

    //- Return the blending field for patchi facei
    virtual scalar patchFacex(const label patchi, const label facei) const
    {
        return x_.boundaryField()[patchi][facei];
    }


public:

    TypeName("cavitating");

    // Constructor
    cavitatingFluidBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    //- Initialize models
    virtual void initializeModels();


    //- Destructor
    virtual ~cavitatingFluidBlastThermo();


    // Member Functions

        //- Correct thermodynamic fields
        virtual void correct();

        //- Add contribution to mixture temperature
        virtual void calculate
        (
            const volScalarField& alpha,
            const volScalarField& he,
            const volScalarField& T,
            volScalarField& Cp,
            volScalarField& Cv,
            volScalarField& mu,
            volScalarField& alphah,
            volScalarField& pXiSum,
            volScalarField& XiSum
        );

        //- Add contribution to speed of sound
        virtual void calculateSpeedOfSound
        (
            const volScalarField& alpha,
            volScalarField& cSqrRhoXiSum
        );

        //- Pre update
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Post update
        virtual void postUpdate();

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;

        //- Return initial energy source
        virtual tmp<volScalarField> initESource() const;


    //- Thermodynamic and transport functions

        //- Enthalpy/Internal energy for a cell [J/kg]
        virtual scalar cellHE
        (
            const scalar T,
            const label celli
        ) const;

        //- Enthalpy/Internal energy for patch and face [J/kg]
        virtual scalar patchFaceHE
        (
            const scalar T,
            const label patchi,
            const label facei
        ) const;

        //- Temperature from enthalpy/internal energy for celli
        virtual scalar cellTHE
        (
            const scalar he,
            const scalar T0,
            const label celli
        ) const;


        using Thermo::updateRho;

        //- Calculate density
        virtual void updateRho(const volScalarField& p);

        //- Calculate density
        virtual void updateRho
        (
            const volScalarField& alpha,
            const volScalarField& p
        );

        //- Calculate pressure for cell [Pa]
        virtual scalar cellpRhoT
        (
            const label celli,
            const bool limit = true
        ) const;

        //- Calculate pressure for cell [Pa]
        virtual scalar patchFacepRhoT
        (
            const label patchi,
            const label facei,
            const bool limit = true
        ) const;

        //- Calculate Gamma = dpde
        virtual tmp<volScalarField> Gamma() const;

        //- Gamma = dpde for cell []
        virtual scalar cellGamma(const label celli) const;

        //- Gamma = dpde for cell []
        virtual scalar patchFaceGamma(const label patchi, const label facei) const;

        //- Derivative of pressure w.r.t density for cell
        virtual scalar celldpdRho(const label celli) const;

        //- Derivative of pressure w.r.t internal energy for cell
        virtual scalar celldpde(const label celli) const;

        //- Derivative of pressure w.r.t temperature for cell
        virtual scalar celldpdT(const label celli) const;

        //- Is the equation of state temperature based
        virtual bool temperatureBased() const
        {
            return Thermo::thermoType1::temperatureBased();
        }

        //- Calculate internal energy from pressure and density
        virtual tmp<volScalarField> calce(const volScalarField& p) const;

        //- Calculate internal energy for a cell
        virtual scalar calcCelle(const scalar p, const label) const;

        //- Calculate the pressure
        virtual tmp<volScalarField> calcp() const;

        //- Calculate the speed of sound
        virtual tmp<volScalarField> calcSpeedOfSound() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cavitatingFluidBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
