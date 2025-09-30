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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeModels
{

/*---------------------------------------------------------------------------*\
                                 Class cavitationModel
\*---------------------------------------------------------------------------*/

class cavitationModel
:
    public phaseChangeModel
{
    // Private data

        //- Is phase 1 the liquid phase
        bool phase1Liquid_;

        //- Saturation vapour pressure
        dimensionedScalar pSat_;

        //- Zero pressure
        dimensionedScalar p0_;


public:

    //- Runtime type information
    TypeName("cavitation");


    // Constructors

        //- Construct for mixture
        cavitationModel
        (
            const word& type,
            const dictionary& dict,
            const volScalarField& alpha1,
            const phaseFluidBlastThermo& thermo1,
            const volScalarField& alpha2,
            const phaseFluidBlastThermo& thermo2
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
            return phase1Liquid_ ? thermo1_.rho() : thermo2_.rho();
        }

        //- Return the internal field of the density of the vapor phase
        const volScalarField::Internal& rhov() const
        {
            return phase1Liquid_ ? thermo2_.rho() : thermo1_.rho();
        }

        //- Return the thermo model of the liquid phase
        const phaseFluidBlastThermo& thermol() const
        {
            return phase1Liquid_ ? thermo1_ : thermo2_;
        }

        //- Return the thermo model of the vapor phase
        const phaseFluidBlastThermo& thermov() const
        {
            return phase1Liquid_ ? thermo2_ : thermo1_;
        }

        //- Return const-access to the saturation vapour pressure
        const dimensionedScalar& pSat() const
        {
            return pSat_;
        }

        //- Return const-access to zero-pressure
        const dimensionedScalar& p0() const
        {
            return p0_;
        }

        //- Return the mass condensation and vaporisation rates
        virtual Pair<tmp<volScalarField::Internal>> mDotCV() const = 0;

        //- Return the mass condensation and vaporisation rates
        virtual Pair<tmp<volScalarField::Internal>> mDots() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseChangeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
