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
    Foam::timeIntegrator

Description
    Base class for time integration of hyperbolic fluxes

SourceFiles
    timeIntegrator.C
    newTimeIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef timeIntegrator_H
#define timeIntegrator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "runTimeSelectionTables.H"
#include "integrationSystem.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class timeIntegrator Declaration
\*---------------------------------------------------------------------------*/

class timeIntegrator
{
protected:
// Protected variables

    //- Reference to mesh
    const fvMesh& mesh_;

    //- Reference to compressible system
    UPtrList<integrationSystem> systems_;


public:

    //- Runtime type information
    TypeName("timeIntegrator");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            timeIntegrator,
            dictionary,
            (const fvMesh& mesh),
            (mesh)
        );

    // Constructor
    timeIntegrator(const fvMesh& mesh);


    //- Destructor
    virtual ~timeIntegrator();


    // Selector
    static autoPtr<timeIntegrator> New(const fvMesh& mesh);


    // Member Functions

        void addSystem(integrationSystem& system);

        virtual void setODEFields(integrationSystem& system) = 0;

        //- Integrate fluxes in time
        virtual void integrate() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
