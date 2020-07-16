/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
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
    Foam::dragModels::Stokes

Description
    Stokes drag model for use in regions with no obvious dispersed phase.

    Reference:
    \verbatim
        "Towards the Numerical Simulation of Multi-scale Two-phase Flows",
        Marschall, H.,
        PhD Thesis, TU München, 2011
    \endverbatim

SourceFiles
    Stokes.C

\*---------------------------------------------------------------------------*/

#ifndef StokesDrag_H
#define StokesDrag_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class Stokes Declaration
\*---------------------------------------------------------------------------*/

class Stokes
:
    public dragModel
{
    // Private data

        //- M coefficient
        const dimensionedScalar dragTime_;


public:

    //- Runtime type information
    TypeName("Stokes");


    // Constructors

        //- Construct from components
        Stokes
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~Stokes();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe(const label, const label) const;

        //- The drag function used in the momentum equation
        virtual tmp<volScalarField> Ki(const label, const label) const;

        //- Drag coefficient
        virtual scalar CdRe
        (
            const label celli,
            const label,
            const label
        ) const;

        //- The drag function used in the momentum equation
        virtual scalar Ki(const label, const label, const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
