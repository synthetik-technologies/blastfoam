/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | Copyright David Blom
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derived work of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "WendlandC4RBFFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBFFunctions
{
    defineTypeNameAndDebug(WendlandC4, 0);
    addToRunTimeSelectionTable(RBFFunction, WendlandC4, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBFFunctions::WendlandC4::WendlandC4(const scalar radius)
:
    RBFFunction(),
    radius_(radius)
{}


Foam::RBFFunctions::WendlandC4::WendlandC4(const dictionary& dict)
:
    RBFFunction(),
    radius_(readScalar(dict.lookup("radius")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBFFunctions::WendlandC4::~WendlandC4()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::RBFFunctions::WendlandC4::evaluate(const scalar dist) const
{
    scalar DByR(dist/radius_);
    return
        dist > radius_
      ? 0.0
      : pow6(1.0 - DByR)*(35.0*sqr(DByR) + 18.0*DByR + 3.0);
}


// ************************************************************************* //
