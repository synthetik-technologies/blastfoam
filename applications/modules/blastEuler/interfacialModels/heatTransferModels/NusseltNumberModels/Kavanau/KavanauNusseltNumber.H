/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::NusseltNumberModels::Kavanau

Description
    Nusselt number of a sphere with the rarefaction correction of Kavanau
    (1955) applied to the continuum correlation of Drake (1961):

        Nu_c = 2 + 0.459 Re^0.55 Pr^0.33
        Nu   = Nu_c/(1 + 3.42 Ma Nu_c/(Re Pr))

    with the particle Mach number Ma = |U_r|/c of the continuous phase.
    Ma/Re = nu/(c d) is evaluated directly, so the correction stays finite
    without slip. The correction was fitted to sphere data for
    0.1 < Ma < 0.69 and 1.75 < Re < 124; it vanishes in the continuum limit
    (Kn = sqrt(pi gamma/2) Ma/Re -> 0).

    References:
    \verbatim
        Kavanau, L.L. (1955). Heat transfer from spheres to a rarefied gas
        in subsonic flow. Transactions of the ASME 77, 617-623.

        Drake, R.M. (1961). Discussion on the paper by G.C. Vliet and
        G. Leppert. Journal of Heat Transfer 83, 170-179.
    \endverbatim

SourceFiles
    KavanauNusseltNumber.C

\*---------------------------------------------------------------------------*/

#ifndef KavanauNusseltNumber_H
#define KavanauNusseltNumber_H

#include "NusseltNumberModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace NusseltNumberModels
{

/*---------------------------------------------------------------------------*\
                          Class Kavanau Declaration
\*---------------------------------------------------------------------------*/

class Kavanau
:
    public NusseltNumberModel
{
public:

    //- Runtime type information
    TypeName("Kavanau");


    // Constructors

        //- Construct from components
        Kavanau
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~Kavanau();


    // Member Functions

        //- Nusselt number from Re, Pr and Ma/Re = nu/(c d)
        static scalar Nu
        (
            const scalar Re,
            const scalar Pr,
            const scalar MaByRe
        );

        //- Nusselt number
        virtual tmp<volScalarField> Nu() const;

        //- Nusselt number for cell
        virtual scalar cellNu(const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace NusseltNumberModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
