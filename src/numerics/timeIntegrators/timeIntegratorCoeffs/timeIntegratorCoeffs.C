/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "timeIntegratorCoeffs.H"
#include "Field.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(timeIntegratorCoeffs, 0);
    defineRunTimeSelectionTable(timeIntegratorCoeffs, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegratorCoeffs::timeIntegratorCoeffs(const label nSteps)
:
    nSteps_(nSteps)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegratorCoeffs::~timeIntegratorCoeffs()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

bool Foam::timeIntegratorCoeffs::update(const label, const scalar)
{
    return false;
}

Foam::List<Foam::label>
Foam::timeIntegratorCoeffs::oldIs(const List<List<scalar>>& as) const
{
    //- Determine if what fields need to be saved
    List<bool> saveOlds(as.size(), false);
    forAll(as, i)
    {
        for (label j = 0; j < as[i].size() - 1; j++)
        {
            saveOlds[j] = saveOlds[j] || mag(as[i][j]) > small;
        }
    }

    label fi = 0;
    List<label> oldIs(as.size(), -1);
    forAll(oldIs, i)
    {
        if (saveOlds[i])
        {
            oldIs[i] = fi++;
        }
    }
    return oldIs;
}


Foam::List<Foam::label>
Foam::timeIntegratorCoeffs::deltaIs(const List<List<scalar>>& bs) const
{
    //- Determine if what fields need to be saved
    List<bool> saveDeltas(bs.size(), false);
    forAll(bs, i)
    {
        for (label j = 0; j < bs[i].size() - 1; j++)
        {
            saveDeltas[j] = saveDeltas[j] || mag(bs[i][j]) > small;
        }
    }

    label fi = 0;
    List<label> deltaIs(bs.size(), -1);
    forAll(deltaIs, i)
    {
        if (saveDeltas[i])
        {
            deltaIs[i] = fi++;
        }
    }
    return deltaIs;
}


void Foam::timeIntegratorCoeffs::setTimeFactors
(
    const List<List<scalar>>& as,
    const List<List<scalar>>& bs,
    List<scalar>& f0,
    List<scalar>& f
) const
{
    f.resize(as.size());
    f0.resize(as.size());
    f0[0] = 0.0;
    f[0] = sum(bs[0]);
    for (label stepi = 1; stepi < as.size(); stepi++)
    {
        scalarList ts(stepi+1, 0.0);
        scalarList dts(stepi, 0.0);
        forAll(dts, i)
        {
            dts[i] = sum(bs[i]);
        }
        ts[1] = dts[0];

        for (label i = 1; i < stepi; i++)
        {
            for (label j = 0; j < as[i].size(); j++)
            {
                ts[i+1] += as[i][j]*ts[j];
            }
            ts[i+1] += dts[i];
        }
        f0[stepi-1] = ts.last() - dts.last();
        f[stepi-1] = f0[stepi-1] + sum(bs[stepi-1]);
    }
}

// ************************************************************************* //
