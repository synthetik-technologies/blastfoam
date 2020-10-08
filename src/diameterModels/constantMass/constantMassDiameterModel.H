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
    Foam::diameterModels::constantMassDiameter

Description
    Constant mass diameter model

SourceFiles
    constantMassDiameter.C

\*---------------------------------------------------------------------------*/

#ifndef constantMassDiameter_H
#define constantMassDiameter_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "diameterModel.H"

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                           Class constantMassDiameter Declaration
\*---------------------------------------------------------------------------*/

class constantMassDiameter
:
    public diameterModel
{
// Private data

    //- Reference density
    dimensionedScalar rho0_;

    //- Reference diameter
    dimensionedScalar d0_;

    //- Reference mass
    dimensionedScalar M0_;

public:

    //- Runtime type information
    TypeName("constantMass");

    // Constructor
    constantMassDiameter
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );


    //- Destructor
    virtual ~constantMassDiameter();

    //- Solve sub-step stepi
    virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
