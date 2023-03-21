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

#include "generalizedPatternSearchMinimizationScheme.H"
#include "SortableList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace minimizationSchemes
{
    defineTypeNameAndDebug(generalizedPatternSearch, 0);
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        generalizedPatternSearch,
        dictionaryMultivariate
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::minimizationSchemes::generalizedPatternSearch::generalizedPatternSearch
(
    const scalarUnivariateEquation& eqns,
    const dictionary& dict
)
:
    minimizationScheme(eqns, dict),
    reflectionCoeff_(dict.lookupOrDefault<scalar>("reflectionCoeff", 1.0)),
    expansionCoeff_(dict.lookupOrDefault<scalar>("expansionCoeff", 2.0)),
    contractionCoeff_(dict.lookupOrDefault<scalar>("contractionCoeff", 0.5))

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::minimizationSchemes::generalizedPatternSearch::minimize
(
    const scalarList& x0,
    const scalarList& xMin,
    const scalarList& xMax,
    const label li
) const
{

}

// ************************************************************************* //
