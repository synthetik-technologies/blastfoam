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
    Foam::detonatingFluidBlastThermo

Description
    Detonating fluid model with options to include activation and afterburn
    models

SourceFiles
    detonatingFluidBlastThermo.C
    detonatingFluidBlastThermos.C

\*---------------------------------------------------------------------------*/

#ifndef detonatingFluidBlastThermo_H
#define detonatingFluidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "activationModel.H"
#include "afterburnModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class detonatingFluidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class detonatingFluidBlastThermo
:
    public Thermo
{
// Private data

    //- Activation model
    autoPtr<activationModel> activation_;

    //- Afterburn model
    autoPtr<afterburnModel> afterburn_;

    //- Correct thermodynamic fields
    void calculate();

    //- Return the blending field for patchi
    virtual tmp<scalarField> x(const label patchi) const
    {
        return activation_->lambdaPow(patchi);
    }

    //- Return the blending field for celli
    virtual scalar cellx(const label celli) const
    {
        return activation_->cellLambdaPow(celli);
    }

    //- Return the blending field for patchi facei
    virtual scalar patchFacex(const label patchi, const label facei) const
    {
        return activation_->patchFaceLambdaPow(patchi, facei);
    }


public:

    TypeName("detonating");

    // Constructor
    detonatingFluidBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    //- Initialize models
    virtual void initializeModels();


    //- Destructor
    virtual ~detonatingFluidBlastThermo();


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

        //- Calculate density
        virtual void updateRho(const volScalarField& p);

        //- Calculate pressure for cell [Pa]
        virtual scalar cellpRhoT
        (
            const label celli,
            const bool limit = true
        ) const;

        //- Calculate Gamma = dpde
        virtual tmp<volScalarField> Gamma() const;

        //- Gamma = dpde for cell []
        virtual scalar cellGamma(const label celli) const;

        //- Derivative of pressure w.r.t density for cell
        virtual scalar celldpdRho(const label celli) const;

        //- Derivative of pressure w.r.t internal energy for cell
        virtual scalar celldpde(const label celli) const;

        //- Calculate internal energy from pressure and density
        virtual tmp<volScalarField> calce(const volScalarField& p) const;

        //- Calculate the pressure
        virtual tmp<volScalarField> calcp() const;

        //- Calculate the speed of sound
        virtual tmp<volScalarField> calcSpeedOfSound() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "detonatingFluidBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
