/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technologies
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
    Foam::scalarLookupTable3D

Description
    Table used to lookup scalar values given a 3D table

SourceFiles
    scalarLookupTable3D.C

\*---------------------------------------------------------------------------*/

#ifndef scalarLookupTable3D_H
#define scalarLookupTable3D_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lookupTable3D.H"
#include "ListOps.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class scalarLookupTable3D Declaration
\*---------------------------------------------------------------------------*/

class scalarLookupTable3D
:
    public lookupTable3D<scalar>
{
// Private data

    //- Include definition of modifying functions
    #include "scalarTableFuncs.H"


public:

    // Constructors

        //- Construct empty
        scalarLookupTable3D();

        //- Construct from dictionary
        scalarLookupTable3D
        (
            const dictionary& dict,
            const word& xName,
            const word& yName,
            const word& zName,
            const word& name
        );

        //- Construct from data
        scalarLookupTable3D
        (
            const Field<scalar>& x,
            const Field<scalar>& y,
            const Field<scalar>& z,
            const Field<Field<Field<scalar>>>& data,
            const word& xMod,
            const word& yMod,
            const word& zMod,
            const word& mod,
            const word& interpolationScheme,
            const bool isReal = true
        );



    //- Destructor
    virtual ~scalarLookupTable3D();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
