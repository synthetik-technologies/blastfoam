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
    Foam::linearAfterburn

Description
    linear afterburn energy

SourceFiles
    linearAfterburn.C

\*---------------------------------------------------------------------------*/

#ifndef linearAfterburn_H
#define linearAfterburn_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "afterburnModel.H"

namespace Foam
{
namespace afterburnModels
{

/*---------------------------------------------------------------------------*\
                           Class linearAfterburn Declaration
\*---------------------------------------------------------------------------*/

class linearAfterburn
:
    public afterburnModel
{
    // Private data

        //- Afterburn energy
        dimensionedScalar Q0_;

        //- Start time
        dimensionedScalar tStart_;

        //- End time
        dimensionedScalar tEnd_;


public:

    //- Runtime type information
    TypeName("linear");

    // Constructor
    linearAfterburn
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );


    //- Destructor
    virtual ~linearAfterburn();


    // Member Functions

        //- Solve sub-step stepi
        virtual void solve
        (
            const label stepi,
            const scalarList& ai,
            const scalarList& bi
        )
        {}

        //- Return energy
        virtual tmp<volScalarField> ESource() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace afterburnModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
