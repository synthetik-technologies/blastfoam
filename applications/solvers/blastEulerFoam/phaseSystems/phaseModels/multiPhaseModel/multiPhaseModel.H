/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2019-04-29 Jeff Heylmun:    Simplified model
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
    Foam::multiPhaseModel

Description
    Base class for a moving phase model. Functions are made so that the class
    can be abstracted to a polydisperse phase.

SourceFiles
    multiPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef multiPhaseModel_H
#define multiPhaseModel_H

#include "fluidPhaseModel.H"
#include "fluxScheme.H"
#include "multiphaseFluidBlastThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class multiPhaseModel Declaration
\*---------------------------------------------------------------------------*/


class multiPhaseModel
:
    public fluidPhaseModel
{
protected:
    // Protected data

        //- Names of components
        wordList phases_;

        //- Primitative variables

            //- Equation of states
            multiphaseFluidBlastThermo& thermo_;

            //- Volume fractions
            PtrList<volScalarField>& alphas_;

            //- Phase densities
            PtrList<volScalarField>& rhos_;


        //- Conserved variables

            //- Phase masses
            PtrList<volScalarField> alphaRhos_;


        //- Fluxes

            //- Volume fraction fluxes
            PtrList<surfaceScalarField> alphaPhis_;

            //- Phase mass fluxes
            PtrList<surfaceScalarField> alphaRhoPhis_;


        //- Calculate new alpha and rho fields
        virtual void calcAlphaAndRho();


public:

    //- Runtime type information
    ClassName("multiphase");


    // Constructor
    multiPhaseModel
    (
        const phaseSystem& fluid,
        const word& phaseName,
        const label index
    );


    //- Destructor
    virtual ~multiPhaseModel();


    // Member Functions

        // integrationSystem functions

            //- Decode primative variables
            virtual void decode();

            //- Encode conserved variables
            virtual void encode();

            //- Update fluxes
            virtual void update();

            //- Solve sub-step stepi
            virtual void solve();

            //- Update viscous terms
            virtual void postUpdate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
