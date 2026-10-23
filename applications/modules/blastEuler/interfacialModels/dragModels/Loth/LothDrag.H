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
    Foam::dragModels::Loth

Description
    Drag of a sphere in compressible and rarefied flow (Loth 2008), as a
    function of the particle Reynolds number Re, the particle Mach number
    Ma = |U_r|/c of the continuous phase and the Knudsen number
    Kn = sqrt(pi gamma/2) Ma/Re:

    Re < 45 (rarefaction dominated):
        Cd = (Cd_Kn,Re + Ma^4 Cd_fm,Re)/(1 + Ma^4)
        Cd_Kn,Re = 24/Re (1 + 0.15 Re^0.687) f_Kn,
        f_Kn = 1/(1 + Kn (2.514 + 0.8 exp(-0.55/Kn)))
        Cd_fm,Re = Cd_fm/(1 + (Cd_fm/1.63 - 1) sqrt(Re/45))
        Cd_fm = (1 + 2 s^2) exp(-s^2)/(s^3 sqrt(pi))
              + (4 s^4 + 4 s^2 - 1) erf(s)/(2 s^4)
              + 2/(3 s) sqrt(pi T_p/T),   s = Ma sqrt(gamma/2)

    Re > 45 (compression dominated):
        Cd = 24/Re (1 + 0.15 Re^0.687) H_M
           + 0.42 C_M/(1 + 42500 G_M/Re^1.16)
        C_M = 5/3 + 2/3 tanh(3 ln(Ma + 0.1))              Ma < 1.45
            = 2.044 + 0.2 exp(-1.8 ln(Ma/1.5)^2)          Ma > 1.45
        G_M = 1 - 1.525 Ma^4                              Ma < 0.89
            = 0.0002 + 0.0008 tanh(12.77 (Ma - 2.02))     Ma > 0.89
        H_M = 1 - 0.258 C_M/(1 + 514 G_M)

    At Ma -> 0 and Kn -> 0 this is the Schiller-Naumann law with the
    Clift-Gauvin term. gamma is the ratio of specific heats of the
    continuous phase (or the optional entry gamma), T_p/T the ratio of the
    dispersed to the continuous phase temperature.

    Reference:
    \verbatim
        Loth, E. (2008). Compressibility and rarefaction effects on drag of
        a spherical particle. AIAA Journal 46(9), 2219-2228.
    \endverbatim

    Usage
    \verbatim
    drag
    (
        (Al in gas)
        {
            type            Loth;
            gamma           1.4;    // optional, else Cp/Cv of the gas
            swarmCorrection
            {
                type        none;
            }
        }
    );
    \endverbatim

SourceFiles
    LothDrag.C

\*---------------------------------------------------------------------------*/

#ifndef LothDrag_H
#define LothDrag_H

#include "dispersedDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                            Class Loth Declaration
\*---------------------------------------------------------------------------*/

class Loth
:
    public dispersedDragModel
{
    // Private data

        //- Ratio of specific heats of the continuous phase (if > 0)
        const scalar gamma_;


    // Private Member Functions

        //- Ratio of specific heats in a cell
        scalar cellGamma(const label celli) const;


public:

    //- Runtime type information
    TypeName("Loth");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Loth
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~Loth();


    // Member Functions

        //- Drag coefficient times Re as a function of Re, Ma, gamma and
        //  the temperature ratio T_p/T
        static scalar CdRe
        (
            const scalar Re,
            const scalar Ma,
            const scalar gamma,
            const scalar TpByT
        );

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;

        //- Drag coefficient
        virtual scalar cellCdRe(const label celli) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
