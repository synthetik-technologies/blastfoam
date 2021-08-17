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
    Foam::dragModels::Chao

Description
    Solid-solid drag model by Chao et. al.

    Reference:
    \verbatim
        Chao, Z., Wang, Y., Jakobsen, J. P., Fernandino, M. & Jaksobsen, H. A.
        (2011).
        "Derivation and validation of a binary multi-fluid Eulerian model for
        fluidized beds."
        Chemical Engineering Science 66(16), 3605-3616
    \endverbatim

SourceFiles
    ChaoDrag.C

\*---------------------------------------------------------------------------*/

#ifndef ChaoDrag_H
#define ChaoDrag_H

#include "dragModel.H"
#include "kineticTheorySystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class Chao Declaration
\*---------------------------------------------------------------------------*/

class Chao
:
    public dragModel
{
    // Private data

        //- Reference to kineticTheorySystem
        const kineticTheorySystem& kineticTheorySystem_;

        //- Coeffient of restitution
        scalar e_;


public:

    //- Runtime type information
    TypeName("Chao");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Chao
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~Chao();


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
