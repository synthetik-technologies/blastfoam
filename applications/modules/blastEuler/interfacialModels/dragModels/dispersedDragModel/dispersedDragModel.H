/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-24 Jeff Heylmun:    Added return functions for acceleration
2025-06-09 Jeff Heylmun:    Added cell based returns
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::dragModels::dispersedDragModel

Description

SourceFiles
    dispersedDragModel.C

\*---------------------------------------------------------------------------*/

#ifndef dispersedDragModel_H
#define dispersedDragModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dragModel.H"
#include "swarmCorrection.H"
#include "phasePair.H"

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                Class dragModels::dispersedDragModel Declaration
\*---------------------------------------------------------------------------*/

class dispersedDragModel
:
    public dragModel
{
protected:

    // Protected data

        //- Swarm correction
        autoPtr<swarmCorrection> swarmCorrection_;


public:

    //- Runtime type information
    TypeName("dispersedDragModel");


    // Constructors

        // Construct with residual constants
        dispersedDragModel
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~dispersedDragModel();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const = 0;

        //- Return the phase-intensive drag coefficient Ki
        //  used in the momentum equations
        //    ddt(alpha1*rho1*U1) + ... = ... alphad*K*(U1-U2)
        //    ddt(alpha2*rho2*U2) + ... = ... alphad*K*(U2-U1)
        virtual tmp<volScalarField> Ki() const;

        //- Return the drag coefficient K
        //  used in the momentum equations
        //    ddt(alpha1*rho1*U1) + ... = ... K*(U1-U2)
        //    ddt(alpha2*rho2*U2) + ... = ... K*(U2-U1)
        virtual tmp<volScalarField> K() const;

        //- Return the drag coefficient Kf
        //  used in the face-momentum equations
        virtual tmp<surfaceScalarField> Kf() const;


        //- Drag coefficient
        virtual scalar cellCdRe(const label) const = 0;

        //- Return the phase-intensive drag coefficient Ki
        //  used in the momentum equations
        //    ddt(alpha1*rho1*U1) + ... = ... alphad*K*(U1-U2)
        //    ddt(alpha2*rho2*U2) + ... = ... alphad*K*(U2-U1)
        virtual scalar cellKi(const label) const;

        //- Return the drag coefficient K
        //  used in the momentum equations
        //    ddt(alpha1*rho1*U1) + ... = ... K*(U1-U2)
        //    ddt(alpha2*rho2*U2) + ... = ... K*(U2-U1)
        virtual scalar cellK(const label) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
