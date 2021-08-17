/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::dragModels::Syamlal

Description
    Solid-sold drag model by Syamlal

    Reference:
    \verbatim
        Syamlal, M. (1987).
        "The particle-particle drag term in a multiparticle model of
        fluidization."
        Topical Report, DOE/MC/21353-2373, NTIS/DE87006500,
        National Technical Information Service, Springfield, VA.

    \endverbatim

SourceFiles
    SyamlalDrag.C

\*---------------------------------------------------------------------------*/

#ifndef SyamlalDrag_H
#define SyamlalDrag_H

#include "dragModel.H"
#include "kineticTheorySystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class Syamlal Declaration
\*---------------------------------------------------------------------------*/

class Syamlal
:
    public dragModel
{
    // Private data

        //- Coefficient of frictional pressure
        dimensionedScalar C1_;

        //- Reference to kineticTheorySystem
        const kineticTheorySystem& kineticTheorySystem_;

        //- Coeffient of restitution
        scalar e_;

        //- Coeffient of friction
        scalar Cf_;


public:

    //- Runtime type information
    TypeName("Syamlal");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Syamlal
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~Syamlal();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe
        (
            const label nodei = -1,
            const label nodej = -1
        ) const;

        //- The drag function used in the momentum equation
        virtual tmp<volScalarField> K
        (
            const label nodei = -1,
            const label nodej = -1
        ) const;

        //- The drag function Kf used in the face-momentum equations
        virtual tmp<surfaceScalarField> Kf
        (
            const label nodei = -1,
            const label nodej = -1
        ) const;

        //- Drag coefficient
        virtual scalar cellCdRe
        (
            const label,
            const label nodei = -1,
            const label nodej = -1
        ) const;

        //- The drag function used in the momentum equation
        virtual scalar cellK
        (
            const label,
            const label nodei = -1,
            const label nodej = -1
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
