/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::embeddedTimeIntegrationSystemBase

Description
    Base class for a solving time dependent ODEs with error calculation

SourceFiles
    embeddedTimeIntegrationSystemBase.C

\*---------------------------------------------------------------------------*/

#ifndef embeddedTimeIntegrationSystemBase_H
#define embeddedTimeIntegrationSystemBase_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegrationSystemBase.H"
#include "embeddedTimeIntegrator.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class embeddedTimeIntegrationSystemBase Declaration
\*---------------------------------------------------------------------------*/

class embeddedTimeIntegrationSystemBase
:
    public timeIntegrationSystemBase
{
protected:
// Protected data


    // Protected member functions

        //- Calculate the error delta given the actual change
        template<template<class> class ListType, class Type>
        Type calcError
        (
            const ListType<Type>& fList
        ) const;

        template<class Type>
        scalar normaliseError
        (
            const Type& y0,
            const Type& y,
            const Type& err,
            const scalar relTol,
            const scalar absTol
        ) const;


public:


    // Constructors

        // Construct from name and object registry
        embeddedTimeIntegrationSystemBase(const word& name);

        // Construct from name and object registry
        embeddedTimeIntegrationSystemBase
        (
            const word& name,
            const objectRegistry& obr
        );

        // Construct from name and timeIntegrator
        embeddedTimeIntegrationSystemBase
        (
            const word& name,
            const timeIntegrator& integrator
        );


    //- Destructor
    virtual ~embeddedTimeIntegrationSystemBase();


    // Member functions

        //- Reset variables when retrying time step
        virtual void reset() = 0;

        //- Return error
        virtual scalar error() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "embeddedTimeIntegrationSystemBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
