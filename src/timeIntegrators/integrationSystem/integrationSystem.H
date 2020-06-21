/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
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
    Foam::integrationSystem

Description
    Base class for a collection of equation of states using a shared pressure
    and velocity (5 equation model)

SourceFiles
    integrationSystem.C

\*---------------------------------------------------------------------------*/

#ifndef integrationSystem_H
#define integrationSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvMesh.H"
#include "Time.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class integrationSystem Declaration
\*---------------------------------------------------------------------------*/

class integrationSystem
:
    public regIOobject
{
protected:
// Protected data

        //- Stored old field indexes
        labelList oldIs_;

        //- Number of stored fields
        label nOld_;


public:


    // Constructor
    integrationSystem(const word& name, const fvMesh& mesh);


    //- Destructor
    virtual ~integrationSystem();


    // Member functions

        //- Decode primative variables
        virtual void decode() = 0;

        //- Encode conserved variables
        virtual void encode() = 0;

        //- Update fluxes
        virtual void update() = 0;

        //- Solve sub-step stepi
        virtual void solve
        (
            const label stepi,
            const scalarList& ai,
            const scalarList& bi
        ) = 0;

        //- Set old lists and fluxes (initialization of fields)
        virtual void setODEFields
        (
            const label nSteps,
            const boolList& storeFields,
            const boolList& storeDeltas
        ) = 0;

        //- Remove stored fields
        virtual void clearODEFields() = 0;


        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
