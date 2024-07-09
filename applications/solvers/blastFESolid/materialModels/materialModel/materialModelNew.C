/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
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

#include "materialModel.H"

// * * * * * * * * * * * * * * Static Data Fucntions * * * * * * * * * * * * //

Foam::autoPtr<Foam::materialModel> Foam::materialModel::New
(
    const dictionary& dict,
    const feMesh1& mesh,
    const pointVectorField& D,
    const pointVectorField& U,
    const bool planeStress,
    const GeoType geoType
)
{
    const word matertialType(dict.lookup("type"));

    if (geoType == LINEAR)
    {
        typename linearConstructorTable::iterator cstrIter =
            linearConstructorTablePtr_->find(matertialType);

        if (cstrIter == linearConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown linear material type " << matertialType << nl
                << "Valid linear material types are:" << nl
                << linearConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }

        return cstrIter()(dict, mesh, D, U, planeStress, geoType);
    }
    else
    {
        typename nonLinearConstructorTable::iterator cstrIter =
            nonLinearConstructorTablePtr_->find(matertialType);

        if (cstrIter == nonLinearConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown non-linear material type " << matertialType << nl
                << "Valid non-linear material types are:" << nl
                << nonLinearConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }

        return cstrIter()(dict, mesh, D, U, planeStress, geoType);
    }
}


// ************************************************************************* //

