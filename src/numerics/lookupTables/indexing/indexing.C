/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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

#include "indexing.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::indexer> Foam::indexer::New
(
    const word& name,
    const List<scalar>& xs
)
{
    if (name == "null")
    {
        return autoPtr<indexer>(new indexers::null(xs));
    }
    else if (name == "uniform")
    {
        return autoPtr<indexer>(new indexers::uniform(xs));
    }
    else if (name =="nonuniform")
    {
        return autoPtr<indexer>(new indexers::nonuniform(xs));
    }
    else
    {
        FatalErrorInFunction
            << "Unknown indexer " << name << nl
            << "Valid indexers are null, uniform, or nonuniform" << endl
            << abort(FatalError);
    }
    return autoPtr<indexer>();
}

Foam::autoPtr<Foam::indexer> Foam::indexer::New
(
    const List<scalar>& xs
)
{
    if (xs.size() < 2)
    {
        return autoPtr<indexer>(new indexers::null(xs));
    }

    bool uniform = true;
    const scalar dx = xs[1] - xs[0];
    for (label i = 2; i < xs.size(); i++)
    {
        if (mag(xs[i] - xs[i-1] - dx) > small)
        {
            uniform = false;
        }
    }

    if (uniform)
    {
        return autoPtr<indexer>(new indexers::uniform(xs));
    }
    else
    {
        return autoPtr<indexer>(new indexers::nonuniform(xs));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::indexers::uniform::findIndex
(
    const scalar x
) const
{
    scalar ij = (x - xs_[0])/(xs_[1] - xs_[0]);
    if (ij <= 0)
    {
        return 0;
    }
    return min(floor(ij), xs_.size() - 2);
}


Foam::label Foam::indexers::nonuniform::findIndex
(
    const scalar x
) const
{
    if (x < xs_[0])
    {
        return 0;
    }

    for (label ij = 0; ij < xs_.size() - 1; ij++)
    {
        if (x > xs_[ij] && x < xs_[ij+1])
        {
            return ij;
        }
    }
    return xs_.size() - 2;
}

// ************************************************************************* //
