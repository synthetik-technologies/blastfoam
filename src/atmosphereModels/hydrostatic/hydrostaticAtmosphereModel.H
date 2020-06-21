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
    Foam::atmosphereModels::hydrostatic

Description
    Initialize atmosphere using a constant density. Pressure is defined so
    the gravitational force and pressure gradient are equal. The reference
    pressure is defined as the pressure at the minimum altitude

SourceFiles
    hydrostaticAtmosphereModel.C

\*---------------------------------------------------------------------------*/

#ifndef hydrostaticAtmosphereModel_H
#define hydrostaticAtmosphereModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "atmosphereModel.H"

namespace Foam
{
namespace atmosphereModels
{

/*---------------------------------------------------------------------------*\
                           Class hydrostatic Declaration
\*---------------------------------------------------------------------------*/

class hydrostatic
:
    public atmosphereModel
{
// Private data

    //- Density value
    dimensionedScalar rhoRef_;

    //- Pressure at reference height
    dimensionedScalar pRef_;


public:

    //- Runtime type information
    TypeName("hydrostatic");

    // Constructor
    hydrostatic(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~hydrostatic();


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
