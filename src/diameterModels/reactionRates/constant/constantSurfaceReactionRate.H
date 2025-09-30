/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
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
    Foam::surfaceReactionRates::constant

Description
    constant reaction rate

SourceFiles
    constantReactionRate.C

\*---------------------------------------------------------------------------*/

#ifndef constantSurfaceReactionRate_H
#define constantSurfaceReactionRate_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceReactionRate.H"

namespace Foam
{
namespace surfaceReactionRates
{

/*---------------------------------------------------------------------------*\
                    Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public surfaceReactionRate
{
// Private data

    //- Rate constant
    dimensionedScalar rate_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from dictionary
        constant(const dictionary& dict);

        //- Construct from runTime and dictionary
        constant(const Time& runTime, const dictionary& dict);

        //- Construct from fvMesh and dictionary
        constant(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~constant();

    //- Construct and return a clone
    virtual autoPtr<surfaceReactionRate> clone() const
    {
        const constant& arr(*this);
        return autoPtr<surfaceReactionRate>(new constant(arr));
    }


    // Public functions

        //- Cell based reaction rate
        virtual scalar k
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const;

        //- Reaction rate
        virtual tmp<volScalarField> k
        (
            const volScalarField& p,
            const volScalarField& T
        ) const;

        //- Reaction rate for a patch
        virtual tmp<scalarField> k
        (
            const fvPatchScalarField& p,
            const fvPatchScalarField& T
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceReactionRates
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
