/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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
    Foam::kineticTheoryModels::viscosityModels::Princeton

Description

    Reference:
    \verbatim
        Agrawal, K., Loezos, P. N., Syamlal, M. & Sundaresan, S. (2001).
        "The role of meso-scale structures in rapid gas–solid flows",
        Journal of Fluid Mechanics, 445, 151-185.
    \endverbatim

SourceFiles
    PrincetonViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef PrincetonViscosity_H
#define PrincetonViscosity_H

#include "kineticTheoryViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class Princeton Declaration
\*---------------------------------------------------------------------------*/

class Princeton
:
    public viscosityModel
{

    typedef
        HashTable
        <
            autoPtr<phasePair>,
            phasePairKey,
            phasePairKey::hash
        >
        phasePairTable;


public:

    //- Runtime type information
    TypeName("Princeton");


    // Constructors

        //- Construct from components
        Princeton
        (
            const dictionary& dict,
            const kineticTheorySystem& kt
        );


    //- Destructor
    virtual ~Princeton();


    // Member functions

        tmp<volScalarField> nu
        (
            const phaseModel& phase,
            const volScalarField& Theta,
            const volScalarField& g0,
            const volScalarField& rho1,
            const volScalarField& da,
            const dimensionedScalar& e
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
