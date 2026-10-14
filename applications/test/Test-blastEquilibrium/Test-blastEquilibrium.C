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

Application
    Test-blastEquilibrium

Description
    Chemical equilibrium at fixed density and internal energy (UV) of a
    mixture of nasa9 species. For each (rho, Tfrozen) pair the absolute
    energy of the initial (frozen) composition at Tfrozen is computed, and the
    equilibrium temperature, ideal-gas pressure and mass fractions at that
    energy and density are printed, with the element and energy errors.

    Test-blastEquilibrium -file thermo.nasa9Gas -species '(N2 O2 NO N O)' \
        -mixture '((N2 0.7669) (O2 0.2331))' -rho '(1 1)' \
        -Tfrozen '(5000 10000)'

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "dictionary.H"
#include "Tuple2.H"
#include "specieBlast.H"
#include "AbelNobel.H"
#include "nasa9BlastThermo.H"
#include "blastChemicalEquilibrium.H"
#include "thermodynamicConstants.H"
#include "cpuTime.H"

using namespace Foam;

typedef nasa9Thermo<AbelNobel<specieBlast>> ThermoType;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("file", "species file");
    argList::addOption("species", "list of species");
    argList::addOption("mixture", "initial composition, (species Y) list");
    argList::addOption("rho", "list of densities [kg/m^3]");
    argList::addOption("Tfrozen", "list of frozen temperatures [K]");
    argList::addOption("coeffs", "equilibrium coefficients dictionary");
    argList args(argc, argv, false, true);
    IOstream::defaultPrecision(17);
    Sout.precision(17);

    IFstream is(args.optionRead<fileName>("file"));
    const dictionary dict(is);
    const hashedWordList species(args.optionRead<wordList>("species"));
    const List<Tuple2<word, scalar>> mixture
    (
        args.optionRead<List<Tuple2<word, scalar>>>("mixture")
    );
    const scalarList rhos(args.optionRead<scalarList>("rho"));
    const scalarList Tfs(args.optionRead<scalarList>("Tfrozen"));
    const dictionary coeffs
    (
        args.optionFound("coeffs")
      ? dictionary(IStringStream(args["coeffs"])())
      : dictionary()
    );

    PtrList<ThermoType> thermos(species.size());
    forAll(species, k)
    {
        thermos.set(k, new ThermoType(dict.subDict(species[k])));
    }

    scalarField Y0(species.size(), 0.0);
    forAll(mixture, i)
    {
        Y0[species[mixture[i].first()]] = mixture[i].second();
    }
    Y0 /= sum(Y0);

    const blastChemicalEquilibrium<ThermoType> eq(species, thermos, coeffs);
    Info<< "elements " << eq.elements() << endl;

    const scalarField b0(eq.elementAbundance(Y0));

    forAll(rhos, s)
    {
        const scalar rho = rhos[s];
        const scalar ea = eq.Ea(Y0, rho, Tfs[s]);

        scalarField Y(Y0);
        scalar T = 3000;
        cpuTime timer;
        const bool ok = eq.solveUV(rho, ea, Y, T);
        const scalar cpu = timer.cpuTimeIncrement();

        scalar p = 0;
        forAll(Y, k)
        {
            if (!eq.condensed()[k])
            {
                p += rho*Y[k]/thermos[k].W()*constant::thermodynamic::RR*T;
            }
        }

        const scalarField b(eq.elementAbundance(Y));
        const scalar bErr = max(mag(b - b0)/max(b0, small));
        const scalar eErr = mag(eq.Ea(Y, rho, T) - ea);

        Info<< "state " << rho << ' ' << ea << ' ' << T << ' ' << p << ' '
            << ok << ' ' << bErr << ' ' << eErr << ' ' << cpu << nl;
        forAll(Y, k)
        {
            Info<< species[k] << ' ' << Y[k] << nl;
        }
    }
    Info<< endl;

    return 0;
}


// ************************************************************************* //
