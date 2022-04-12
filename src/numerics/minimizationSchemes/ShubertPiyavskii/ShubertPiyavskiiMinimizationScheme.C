/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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

#include "ShubertPiyavskiiMinimizationScheme.H"
#include "ListOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ShubertPiyavskiiMinimizationScheme, 0);
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        ShubertPiyavskiiMinimizationScheme,
        dictionaryUnivariate
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        ShubertPiyavskiiMinimizationScheme,
        dictionaryZero
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        ShubertPiyavskiiMinimizationScheme,
        dictionaryOne
    );
    addToRunTimeSelectionTable
    (
        univariateMinimizationScheme,
        ShubertPiyavskiiMinimizationScheme,
        dictionaryTwo
    );
}

Foam::vector2D Foam::ShubertPiyavskiiMinimizationScheme::intersection
(
    const vector2D& A,
    const vector2D& B
) const
{
    scalar t = (A.y() - B.y() - l_*(A.x() - B.x()))/(2.0*l_);
    return vector2D(A.x() + t, A.y() - t*l_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ShubertPiyavskiiMinimizationScheme::ShubertPiyavskiiMinimizationScheme
(
    const scalarUnivariateEquation& eqn,
    const dictionary& dict
)
:
    univariateMinimizationScheme(eqn, dict),
    l_
    (
        dict.lookupOrDefault<scalar>
        (
            "l",
            (eqn_.upper() - eqn_.lower())/2.0
        )
    )
{
    checkY_ = true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::ShubertPiyavskiiMinimizationScheme::minimize
(
    const scalar x,
    const scalar x1,
    const scalar x2,
    const label li
) const
{
    vector2D A(x1, eqn_.fx(x1, li));
    vector2D B(x2, eqn_.fx(x2, li));
    vector2D P(x, eqn_.fx(x, li));
    vector2D P_prev(intersection(A, P));
    vector2D P_next(intersection(P, B));
    DynamicList<scalar> ptsx({A.x(), P_prev.x(), P.x(), P_next.x(), B.x()});
    DynamicList<scalar> ptsy({A.y(), P_prev.y(), P.y(), P_next.y(), B.y()});

    scalar delta = great;
    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        label i = findMin(ptsy);
        P = vector2D(ptsx[i], eqn_.fx(ptsx[i], li));
        delta = mag(P.y() - ptsy[i]);

        P_prev = vector2D(ptsx[i-1], ptsy[i-1]);
        P_next = vector2D(ptsx[i+1], ptsy[i+1]);

        P_prev = intersection(P_prev, P);
        P_next = intersection(P, P_next);

        if
        (
            convergedX(P_next.x(), P_prev.x())
         && convergedY(P_next.y(), P_prev.y())
        )
        {
            break;
        }

        List<scalar> ptsxTmp = ptsx;
        List<scalar> ptsyTmp = ptsy;
        ptsx.resize(ptsx.size() + 2);
        ptsy.resize(ptsx.size());

        ptsx[i] = P_prev.x();
        ptsx[i+1] = P.x();
        ptsx[i+2] = P_next.x();

        ptsy[i] = P_prev.y();
        ptsy[i+1] = P.y();
        ptsy[i+2] = P_next.y();
        for (label j = i+1; j < ptsxTmp.size(); j++)
        {
            ptsx[j+2] = ptsxTmp[j];
            ptsy[j+2] = ptsyTmp[j];
        }
        printStepInformation(P.x());
    }

    return printFinalInformation(P.x());
}

// ************************************************************************* //
