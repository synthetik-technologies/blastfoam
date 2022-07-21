/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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

#include "lookupTables1D.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
Foam::scalar
Foam::lookupTable1D<Foam::scalar>::reverseLookup(const scalar& fin) const
{
#ifdef FULL_DEBUG
    if (!mod_.valid())
    {
        FatalErrorInFunction
            << "Try to interpolate data that has not been set."
            << abort(FatalError);
    }
#endif

    scalar f(mod_()(fin));
    if (f < data_[0])
    {
        index_ = 0;
    }
    if (f > data_.last())
    {
        index_ = data_.size() - 2;
    }
    else
    {
        for (index_ = 0; index_ < data_.size(); index_++)
        {
            if (f < data_[index_])
            {
                index_--;
                break;
            }
        }
    }

    const scalar& fm(data_[index_]);
    const scalar& fp(data_[index_+1]);

    scalar w = interpolationWeight1D::linearWeight(f, fm, fp);

    return modX_->inv
    (
        xModValues_[index_]
      + w*(xModValues_[index_+1] - xModValues_[index_])
    );
}

// ************************************************************************* //
