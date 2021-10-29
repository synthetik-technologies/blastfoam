/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
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

#include "immersedShape.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
Foam::autoPtr<Foam::immersedShape> Foam::immersedShape::New
(
    const polyMesh& mesh,
    const immersedBoundaryObject& ibo,
    const dictionary& dict
)
{
    const word shapeType(dict.lookup("shape"));
    Info<< "Selecting shape: " << shapeType << endl;

    if (mesh.nGeometricD() == 3)
    {
        threeDConstructorTable::iterator cstrIter =
            threeDConstructorTablePtr_->find(shapeType);

        if (cstrIter == threeDConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown 3-D immersedShape type " << shapeType
                << endl << endl
                << "Valid 3-D immersedShape types : " << endl
                << threeDConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<immersedShape>(cstrIter()(mesh, ibo, dict));
    }
    else if (mesh.nGeometricD() == 2)
    {
        twoDConstructorTable::iterator cstrIter =
            twoDConstructorTablePtr_->find(shapeType);

        if (cstrIter == twoDConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown 2-D immersedShape type " << shapeType
                << endl << endl
                << "Valid 2-D immersedShape types : " << endl
                << twoDConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<immersedShape>(cstrIter()(mesh, ibo, dict));
    }
    else
    {
        FatalErrorInFunction
            << "1-D meshes are not valid with immersed objects"
            << endl;
    }
    return autoPtr<immersedShape>();
}


// ************************************************************************* //
