/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::blendingMethods::noBlending

Description
    No blending: in each pair, the model for the phase that is dispersed is
    used. Either the dispersed phase (dispersedPhase) or the continuous
    phase is named. The continuous phase can be one phase (continuousPhase)
    or a list in order of precedence (continuousPhases): in a pair, the
    phase that comes first in the list is the continuous one, e.g.

        continuousPhases (gas compBAl);

    makes compBAl continuous in the pair (Al, compBAl) and dispersed in
    (compBAl, gas). A pair with neither phase in the list uses only its
    unordered model (as with a single continuousPhase).

SourceFiles
    noBlending.C

\*---------------------------------------------------------------------------*/

#ifndef noBlending_H
#define noBlending_H

#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{

/*---------------------------------------------------------------------------*\
                         Class noBlending Declaration
\*---------------------------------------------------------------------------*/

class noBlending
:
    public blendingMethod
{
    // Private Data

        //- Names of the continuous phases in order of precedence
        wordList continuousPhases_;

        //- Name of the dispersed phase
        const word dispersedPhase_;


    // Private Member Functions

        //- Precedence of a phase as the continuous phase (lower first;
        //  labelMax if not listed)
        label rank(const word& phaseName) const;


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from a dictionary and two phases
        noBlending
        (
            const dictionary& dict,
            const wordList& phaseNames
        );


    //- Destructor
    ~noBlending();


    // Member Functions

        //- Factor for first phase
        virtual tmp<volScalarField> f1
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Factor for first phase for a cell
        virtual scalar f1i
        (
            const label celli,
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Factor for second phase
        virtual tmp<volScalarField> f2
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Factor for second phase for a cell
        virtual scalar f2i
        (
            const label celli,
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
