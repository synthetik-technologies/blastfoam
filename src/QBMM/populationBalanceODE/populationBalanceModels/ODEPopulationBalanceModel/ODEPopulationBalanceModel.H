/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 Alberto Passalacqua
     \\/     M anipulation  |
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
    Foam::ODEPopulationBalanceModel

Description
    Abstract class for the run-time selection of the population balance model.

SourceFiles
    ODEPopulationBalanceModel.C

\*---------------------------------------------------------------------------*/

#ifndef ODEPopulationBalanceModel_H
#define ODEPopulationBalanceModel_H

#include "IOdictionary.H"
#include "fvCFD.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "integrationSystem.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ODEPopulationBalanceModel Declaration
\*---------------------------------------------------------------------------*/

class ODEPopulationBalanceModel
:
    public integrationSystem
{
    // Private data

        //- Name of the ODEPopulationBalanceModel
        const word name_;


protected:

    // Protected data

        //- Reference to populationBalanceProperties
        //  Non-const so ode solver can read on modification
        dictionary& populationBalanceProperties_;

        //- Fluid-phase face-normal velocity
        const surfaceScalarField& phi_;


public:

    //- Runtime type information
    TypeName("ODEPopulationBalanceModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ODEPopulationBalanceModel,
        dictionary,
        (
            const word& name,
            const dictionary& dict,
            const surfaceScalarField& phi
        ),
        (name, dict, phi)
    );


    // Constructors

        //- Construct from components
        ODEPopulationBalanceModel
        (
            const word& name,
            const dictionary& dict,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construct
        ODEPopulationBalanceModel(const ODEPopulationBalanceModel&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<ODEPopulationBalanceModel> New
        (
            const word& name,
            const dictionary& dict,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~ODEPopulationBalanceModel();


    // Member Functions

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo() const = 0;

        //- Return the courant number based on abscissae
        virtual scalar CoNum() const = 0;

        //- Update fields
        virtual void update()
        {
            return;
        }

        //- Solve sub-step stepi
        virtual void solve() = 0;

        //- Set old lists and fluxes (initialization of fields)
        virtual void postUpdate() = 0;

        //- Remove stored fields
        virtual void clearODEFields() = 0;


    // regIOobject functions

        //- Read
        virtual bool readIfModified()
        {
            return true;
        }

        //- Read
        virtual bool writeData(Ostream&) const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ODEPopulationBalanceModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
