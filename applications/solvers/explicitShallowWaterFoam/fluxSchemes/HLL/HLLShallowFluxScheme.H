/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2019-10-21  Jeff Heylmun:   Moved from rhoCentralFoam to runtime selectable
                            method.
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
    Foam::shallowFluxSchemes::HLL

Description

SourceFiles
    HLLShallowFluxScheme.C

\*---------------------------------------------------------------------------*/

#ifndef HLLShallowFluxScheme_H
#define HLLShallowFluxScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "shallowFluxScheme.H"

namespace Foam
{
namespace shallowFluxSchemes
{

/*---------------------------------------------------------------------------*\
                           Class HLL Declaration
\*---------------------------------------------------------------------------*/

class HLL
:
    public shallowFluxScheme
{
protected:
    // Private functions

        //- Calcualte fluxes
        virtual void calculateFluxes
        (
            const scalar& hOwn, const scalar& hNei,
            const vector& UOwn, const vector& UNei,
            const vector& Sf,
            scalar& phi,
            scalar& hPhi,
            vector& hUPhi,
            const label facei, const label patchi = -1
        );


public:

    //- Runtime type information
    TypeName("HLL");

    // Constructor
    HLL
    (
        surfaceScalarField& phi,
        surfaceScalarField& hPhi,
        surfaceVectorField& hUPhi,
        const dimensionedVector& g
    );


    //- Destructor
    virtual ~HLL();


    // Member Functions

        //- Clear savedFields
        virtual void clear();

        //- Allocate saved fields
        virtual void createSavedFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace shallowFluxSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
