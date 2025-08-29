/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::phaseChangeModels::cavitationModel

Description
    Abstract base class for cavitation models

SourceFiles
    cavitationModel.C

\*---------------------------------------------------------------------------*/

#ifndef cavitationModel_H
#define cavitationModel_H

#include "phaseChangeModel.H"
#include "saturationPressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                                 Class cavitationModel
\*---------------------------------------------------------------------------*/

class cavitationModel
:
    public phaseChangeModel
{
public:

    //- Zero pressure
    static const dimensionedScalar zeroP;


private:
    // Private data

        //- Is phase 1 the liquid phase
        bool phase1Liquid_;

        //- Phase 1 density
        const volScalarField& rho1_;

        //- Phase 2 density
        const volScalarField& rho2_;

        //- Phase 1 temperature
        const volScalarField& T1_;

        //- Phase 2 temperature
        const volScalarField& T2_;

        //- Saturation pressure model
        autoPtr<saturationPressureModel> pSat_;


public:

    //- Runtime type information
    TypeName("cavitationModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            cavitationModel,
            phaseThermo,
            (
                const dictionary& dict,
                const volScalarField& alpha1,
                const phaseFluidBlastThermo& thermo1,
                const volScalarField& alpha2,
                const phaseFluidBlastThermo& thermo2
            ),
            (dict, alpha1, thermo1, alpha2, thermo2)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            cavitationModel,
            thermo,
            (
                const dictionary& dict,
                const volScalarField& alpha1,
                const fluidThermo& thermo1,
                const volScalarField& alpha2,
                const fluidThermo& thermo2
            ),
            (dict, alpha1, thermo1, alpha2, thermo2)
        );

    // Constructors

        //- Construct for mixture
        cavitationModel
        (
            const word& type,
            const dictionary& dict,
            const volScalarField& alpha1,
            const volScalarField& rho1,
            const volScalarField& T1,
            const volScalarField& alpha2,
            const volScalarField& rho2,
            const volScalarField& T2
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<cavitationModel> New
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const phaseFluidBlastThermo& thermo1,
            const volScalarField& alpha2,
            const phaseFluidBlastThermo& thermo2
        );

        //- Return a reference to the selected phaseChange model
        static autoPtr<cavitationModel> New
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const fluidThermo& thermo1,
            const volScalarField& alpha2,
            const fluidThermo& thermo2
        );

    //- Destructor
    virtual ~cavitationModel()
    {}


    // Member Functions

        //- Return liquid phase volume fraction
        const volScalarField::Internal& alphal() const
        {
            return phase1Liquid_ ? alpha1_ : alpha2_;
        }

        //- Return vapor volume fraction
        const volScalarField::Internal& alphav() const
        {
            return phase1Liquid_ ? alpha2_ : alpha1_;
        }

        //- Return the internal field of the density of the liquid phase
        const volScalarField::Internal& rhol() const
        {
            return phase1Liquid_ ? rho1_ : rho2_;
        }

        //- Return the internal field of the density of the vapor phase
        const volScalarField::Internal& rhov() const
        {
            return phase1Liquid_ ? rho2_ : rho2_;
        }

        //- Return the phase1 saturation pressure
        tmp<volScalarField::Internal> pSat1() const
        {
            return pSat_->pSat(T1_);
        }

        //- Return phase2 saturation pressure
        tmp<volScalarField::Internal> pSat2() const
        {
            return pSat_->pSat(T2_);
        }

        //- Return the saturation liquid pressure
        tmp<volScalarField::Internal> pSatl() const
        {
            return phase1Liquid_ ? pSat1() : pSat2();
        }

        //- Return the saturation vapour pressure
        tmp<volScalarField::Internal> pSatv() const
        {
            return phase1Liquid_ ? pSat2() : pSat1();
        }


        //- Return the mass condensation and vaporisation rates
        virtual Pair<tmp<volScalarField::Internal>> mDotCV() const = 0;

        //- Return the mass condensation and vaporisation rates
        virtual Pair<tmp<volScalarField::Internal>> mDots() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
