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
    Foam::blastChemicalEquilibrium

Description
    Chemical equilibrium of a multicomponent gas at fixed density and
    absolute specific internal energy (UV problem), or at fixed density and
    temperature (TV problem), by the element-potential method.

    Gas species are ideal in the chemical potentials,
        mu_k/(R_u T) = g_k + ln(c_k R_u T/p_std),   g_k = W_k G_std,k/(R_u T),
    with c_k the molar concentration. Condensed species (name with a phase
    suffix in parentheses, e.g. C(gr), Al2O3(s)) are pure phases,
    mu_c/(R_u T) = g_c. At fixed T and V the gas amounts follow from the
    element potentials lambda,
        n_k = V p_std/(R_u T) exp(-g_k + sum_i a_ik lambda_i),
    and lambda maximises the concave dual
        b.lambda - sum_k n_k(lambda)
    subject to sum_i a_ic lambda_i <= g_c for each condensed species, whose
    Lagrange multipliers are the condensed amounts. This is solved by Newton
    iterations on the KKT system with an active set of condensed species and
    a backtracking line search on the Lagrangian; the KKT system is scaled
    by the element abundances and solved with a pseudo-inverse, so that
    directions not fixed by significant species (elements held in a
    condensed phase) get no step. The temperature at the given energy is
    found by a safeguarded Newton iteration on the equilibrium energy.

    The energy of each species is that of the thermo model, Es + Hf (for the
    Abel-Nobel gas, h - R T; the covolume enters the pressure only), so that
    the equilibrium composition conserves the energy of the flow solver.

    Elements are parsed from the species names (e.g. Al2O3(s) = 2 Al + 3 O);
    the element composition of a species may be given explicitly with
    \verbatim
        elements
        {
            <specie> ((<element> <count>) ...);
        }
    \endverbatim
    in the coefficient dictionary.

SourceFiles
    blastChemicalEquilibrium.C

\*---------------------------------------------------------------------------*/

#ifndef blastChemicalEquilibrium_H
#define blastChemicalEquilibrium_H

#include "speciesTable.H"
#include "PtrList.H"
#include "scalarField.H"
#include "boolList.H"
#include "labelList.H"
#include "dictionary.H"
#include "scalarMatrices.H"
#include "Tuple2.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class blastChemicalEquilibrium Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class blastChemicalEquilibrium
{
    // Private Data

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& thermos_;

        //- Number of species
        const label nSpecie_;

        //- Element names
        wordList elements_;

        //- Number of atoms of element i in specie k, a_[k][i]
        List<scalarList> a_;

        //- Is the specie a condensed (pure) phase
        boolList condensed_;

        //- Relative tolerance (element balance and energy)
        scalar tolerance_;

        //- Maximum number of iterations (each of the T and the lambda loops)
        label maxIter_;

        //- Maximum change of an element potential in a Newton step
        scalar maxStep_;

        //- Temperature bounds of the energy iteration
        scalar Tmin_;
        scalar Tmax_;

        //- Print the iterations
        bool debug_;


    // Private Member Functions

        //- Element composition parsed from a specie name
        static List<Tuple2<word, scalar>> parseFormula(const word& name);

        //- Dimensionless standard Gibbs energy g_k = W_k G_std,k/(R_u T)
        void gibbs(const scalar T, scalarField& g) const;

        //- Gas amounts [kmol/kg] for the element potentials
        void gasAmounts
        (
            const scalar lnc0,
            const scalarField& g,
            const boolList& active,
            const scalarField& lambda,
            scalarField& n
        ) const;

        //- Initial element potentials from a composition (least squares)
        void initialLambda
        (
            const scalar T,
            const scalar V,
            const scalarField& n0,
            const scalarField& b,
            const boolList& activeEl,
            const boolList& active,
            scalarField& lambda
        ) const;

        //- Solve the KKT system for the given right-hand side
        //  (element rows r, condensed rows rc); returns false if singular
        bool solveKKT
        (
            const scalarField& n,
            const scalarField& b,
            const boolList& activeEl,
            const labelList& C,
            const scalarField& r,
            const scalarField& rc,
            scalarField& dLambda,
            scalarField& nC
        ) const;

        //- TV equilibrium for the element abundances b [kmol/kg];
        //  lambda is the initial guess on input. Returns n [kmol/kg],
        //  dn/dT and whether converged
        bool solveTV
        (
            const scalar T,
            const scalar V,
            const scalarField& b,
            const boolList& activeEl,
            const boolList& active,
            scalarField& lambda,
            DynamicList<label>& C,
            scalarField& n,
            scalarField& dndT
        ) const;

        //- Active elements and species of the element abundances b
        void activeSets
        (
            const scalarField& b,
            boolList& activeEl,
            boolList& active
        ) const;


public:

    // Constructors

        //- Construct from the species, their thermo and the coefficients
        blastChemicalEquilibrium
        (
            const speciesTable& species,
            const PtrList<ThermoType>& thermos,
            const dictionary& dict
        );


    // Member Functions

        //- Element names
        const wordList& elements() const
        {
            return elements_;
        }

        //- Atoms of each element per specie
        const List<scalarList>& atoms() const
        {
            return a_;
        }

        //- Is the specie condensed
        const boolList& condensed() const
        {
            return condensed_;
        }

        //- Element abundances [kmol/kg] of the mass fractions Y
        tmp<scalarField> elementAbundance(const scalarField& Y) const;

        //- Absolute specific internal energy [J/kg] at (rho, T, Y)
        scalar Ea(const scalarField& Y, const scalar rho, const scalar T)
            const;

        //- Equilibrium mass fractions at fixed rho and T; Y is the initial
        //  guess on input. Returns whether converged
        bool solveTV(const scalar rho, const scalar T, scalarField& Y) const;

        //- Equilibrium mass fractions and temperature at fixed rho and
        //  absolute energy ea; Y and T are the initial guesses on input.
        //  Returns whether converged
        bool solveUV
        (
            const scalar rho,
            const scalar ea,
            scalarField& Y,
            scalar& T
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "blastChemicalEquilibrium.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
