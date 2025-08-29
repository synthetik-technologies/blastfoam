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
    Foam::phaseChangeModel

Description

SourceFiles
    phaseChangeModel.C
    phaseChangeModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeModel_H
#define phaseChangeModel_H

#include "fluidBlastThermo.H"
#include "phaseFluidBlastThermo.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class phaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeModel
{
protected:

    // Protected data

        //- Constructing dictionary
        const dictionary dict_;

        //- First phase volume faction
        const volScalarField& alpha1_;

        //- Second phase volume faction
        const volScalarField& alpha2_;


public:

    //- Runtime type information
    TypeName("phaseChangeModel");

    //- Default phaseChangeProperties dictionary name
    static const word phaseChangePropertiesName;


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseChangeModel,
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
            phaseChangeModel,
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

        //- Construct from dictionary
        phaseChangeModel
        (
            const word& type,
            const dictionary& dict,
            const volScalarField& alpha1,
            const volScalarField& alpha2
        );

        //- Disallow default bitwise copy construction
        phaseChangeModel(const phaseChangeModel&);


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<phaseChangeModel> New
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const phaseFluidBlastThermo& thermo1,
            const volScalarField& alpha2,
            const phaseFluidBlastThermo& thermo2
        );

        //- Return a reference to the selected phaseChange model
        static autoPtr<phaseChangeModel> New
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const fluidThermo& thermo1,
            const volScalarField& alpha2,
            const fluidThermo& thermo2
        );


    //- Destructor
    virtual ~phaseChangeModel()
    {}


    // Member Functions

        //- Return the mass condensation and vaporisation rates
        virtual Pair<tmp<volScalarField::Internal>> mDots() const = 0;

        //- Correct the phaseChange model
        virtual void correct() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseChangeModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
