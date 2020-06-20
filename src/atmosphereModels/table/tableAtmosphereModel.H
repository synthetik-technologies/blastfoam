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
    Foam::atmosphereModels::table

Description
    Initialize atmosphere using pressure and temperature tables based on
    height. An ideal gas is assumed.

SourceFiles
    tableAtmosphereModel.C

\*---------------------------------------------------------------------------*/

#ifndef tableAtmosphereModel_H
#define tableAtmosphereModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lookupTable1D.H"
#include "atmosphereModel.H"

namespace Foam
{
namespace atmosphereModels
{

/*---------------------------------------------------------------------------*\
                           Class table Declaration
\*---------------------------------------------------------------------------*/

class table
:
    public atmosphereModel
{
// Private data

    //- Pressure table
    lookupTable1D pTable_;

    //- Temperature table
    lookupTable1D TTable_;

    //- Molecular weight
    scalar W_;


public:

    //- Runtime type information
    TypeName("table");

    // Constructor
    table(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~table();


    // Member Functions

        //- Create atmosphere
        virtual void createAtmosphere
        (
            volScalarField& p,
            volScalarField& rho
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace atmosphereModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
