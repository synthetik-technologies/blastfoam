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

#include "pointFiniteElement.H"
#include "addToRunTimeSelectionTable.H"
#include "addToRunTimeSelectionMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{
    defineTypeNameAndDebug(point, 0);
    addToRunTimeSelectionTable(finiteElement, point, type);
    addNamedToRunTimeSelectionTable(finiteElement, point, type, pt);

    addToRunTimeSelectionTable(finiteElement, point, typeOrder);
    addToRunTimeSelectionMap(finiteElement, point, msh, 0);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::finiteElements::point::point(const label order)
:
    FiniteElement<ElementType::PT>(order, 1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::finiteElements::point::~point()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::finiteElements::point::getNodes
(
    const List<vector>& verts
) const
{
    return tmp<pointField>(new pointField(verts));
}


Foam::scalarList Foam::finiteElements::point::calcShape
(
    const vector& pt
) const
{
    return scalarList({1.0});
}


Foam::scalarRectangularMatrix Foam::finiteElements::point::calcDShape
(
    const vector& pt
) const
{
    return scalarRectangularMatrix(1, 1, 0.0);
}


void Foam::finiteElements::point::vtkData
(
    labelList& data,
    label& start,
    const labelList& labels
) const
{
    data[start++] = labels[0];
}

// ************************************************************************* //

