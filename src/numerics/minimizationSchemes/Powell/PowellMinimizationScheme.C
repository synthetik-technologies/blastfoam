/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
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

\*---------------------------------------------------------------------------*/

#include "PowellMinimizationScheme.H"
#include "SortableList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace minimizationSchemes
{
    defineTypeNameAndDebug(Powell, 0);
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        Powell,
        dictionaryMultivariate
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::minimizationSchemes::Powell::Powell
(
    const scalarUnivariateEquation& eqns,
    const dictionary& dict
)
:
    basis(eqns, dict),
    dirs_(eqns.nVar(), scalarList(eqns.nVar(), 0.0))
{}


Foam::minimizationSchemes::Powell::Powell
(
    const scalarUnivariateEquation& eqns,
    const Powell& solver
)
:
    basis(eqns, solver),
    dirs_(eqns.nVar(), scalarList(eqns.nVar(), 0.0))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::minimizationSchemes::Powell::minimize
(
    const scalarList& x0,
    const scalarList& xMin,
    const scalarList& xMax,
    const label li
) const
{
    tmp<scalarField> txNew(new scalarField(x0));
    scalarField& xNew = txNew.ref();
    scalarField xOld(xNew);
    scalarField delta(xNew);

    forAll(dirs_, i)
    {
        dirs_[i] = 0.0;
        dirs_[i][i] = 1.0;
    }

    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        xOld = xNew;
        forAll(xNew, diri)
        {
            if (norm(dirs_[diri]) > small)
            {
                eqns_.limit(xNew);
                lineSearcher().searchDir(xNew, dirs_[diri], li, xNew);
            }
        }
        for (label diri = 0; diri < dirs_.size()-1; diri++)
        {
            dirs_[diri] = dirs_[diri+1];
        }
        delta = xNew - xOld;
        dirs_.last() = delta;

        if (norm(delta) > small)
        {
            eqns_.limit(xNew);
            lineSearcher().searchDir(xNew, delta, li, xNew);
            delta = xNew - xOld;
        }

        if (convergedXScale(delta, xNew))
        {
            break;
        }
        printStepInformation(xNew);
    }
    printFinalInformation(xNew);
    return txNew;
}

// ************************************************************************* //
