/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technologies
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
    Foam::diameterModel

Description
    Base model to describe the diamter of a phase

SourceFiles
    diameterModel.C
    newDiameterModel.C

\*---------------------------------------------------------------------------*/

#ifndef diameterModel_H
#define diameterModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "integrationSystem.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class diameterModel Declaration
\*---------------------------------------------------------------------------*/

class diameterModel
:
    public integrationSystem
{
protected:

    // Protected data

        //- Diameter
        volScalarField d_;

        //- Does the diameter field need to be read
        virtual bool requireRead()
        {
            return false;
        }


public:

    //- Runtime type information
    ClassName("diameterModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            diameterModel,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const word& phaseName
            ),
            (mesh, dict, phaseName)
        );

    // Constructor
    diameterModel
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );


    //- Destructor
    virtual ~diameterModel();


    // Selectors

        static autoPtr<diameterModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName
        );


    // Member Functions

        //- Update
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Solve sub-step stepi given pressure and temperature fields
        virtual void solve(const volScalarField& p, const volScalarField& T);

        //- Post update
        virtual void postUpdate();

        //- Remove stored fields
        virtual void clearODEFields();

        //- Return if the model requires inputs
        virtual bool requireInput() const
        {
            return false;
        }

        //- Return effective diameter
        const volScalarField& d() const
        {
            return d_;
        }

        //- Return surface area
        virtual tmp<volScalarField> A() const;

        //- Return volume
        virtual tmp<volScalarField> V() const;

        //- Return change in radius per unit time
        //  Only used in mass transfer
        virtual tmp<volScalarField> dMdt() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
