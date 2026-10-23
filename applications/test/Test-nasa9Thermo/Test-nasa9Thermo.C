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
    Test-nasa9Thermo

Description
    Prints Cp/R, H/R and S/R (standard state, per kmol) of the species of a
    nasa9 species file, and of a mixture given as mass fractions, at a list
    of temperatures, for comparison with the source data.

    Test-nasa9Thermo -file thermo.nasa9Gas -species '(N2 O2 O)' \
        -mixture '((N2 0.7669) (O2 0.2331))' -T '(200 300 1000 6000 20000)'

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "dictionary.H"
#include "Tuple2.H"
#include "specieBlast.H"
#include "AbelNobel.H"
#include "nasa9BlastThermo.H"

using namespace Foam;

typedef nasa9Thermo<AbelNobel<specieBlast>> ThermoType;

void print(const word& name, const ThermoType& t, const List<scalar>& Ts)
{
    // Per kmol: R = Ru/W
    const scalar R = t.R();
    forAll(Ts, i)
    {
        const scalar T = Ts[i];
        const scalar Ha = t.Hat(T);
        const scalar S = (Ha - t.Gstd(T))/T;
        Info<< name << ' ' << T << ' ' << t.Cpt(T)/R << ' ' << Ha/R << ' '
            << S/R << endl;
    }
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("file", "species file");
    argList::addOption("species", "list of species to print");
    argList::addOption("mixture", "list of (species massFraction)");
    argList::addOption("T", "list of temperatures [K]");
    argList args(argc, argv, false, true);

    IFstream is(args.optionRead<fileName>("file"));
    const dictionary dict(is);
    const List<scalar> Ts(args.optionRead<List<scalar>>("T"));

    Info<< "# name T Cp/R H/R S/R" << endl;

    if (args.optionFound("species"))
    {
        const wordList species(args.optionRead<wordList>("species"));
        forAll(species, i)
        {
            const ThermoType t(dict.subDict(species[i]));
            print(species[i], t, Ts);
        }
    }

    if (args.optionFound("mixture"))
    {
        const List<Tuple2<word, scalar>> Y
        (
            args.optionRead<List<Tuple2<word, scalar>>>("mixture")
        );
        autoPtr<ThermoType> mix;
        forAll(Y, i)
        {
            const ThermoType t(dict.subDict(Y[i].first()));
            if (mix.valid())
            {
                mix() += Y[i].second()*t;
            }
            else
            {
                mix.set(new ThermoType("mixture", Y[i].second()*t));
            }
        }
        print("mixture", mix(), Ts);
    }

    return 0;
}


// ************************************************************************* //
