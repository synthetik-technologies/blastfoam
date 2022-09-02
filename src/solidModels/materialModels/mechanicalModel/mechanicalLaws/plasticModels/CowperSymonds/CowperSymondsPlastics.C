/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
31-03-2022 Synthetik Applied Technologies: Added CowperSymondsPlastic
-------------------------------------------------------------------------------
License
    This file is a derivative work of foam-extend.

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

#include "CowperSymondsPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "linearPlasticModel.H"
#include "nonLinearPlasticModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    typedef CowperSymondsPlastic<linearPlasticModel> linearCowperSymondsPlastic;
    addNamedToRunTimeSelectionTable
    (
        mechanicalLaw,
        linearCowperSymondsPlastic,
        linGeomMechLaw,
        linearCowperSymondsPlastic
    );
    typedef CowperSymondsPlastic<nonLinearPlasticModel> nonLinearCowperSymondsPlastic;
    addNamedToRunTimeSelectionTable
    (
        mechanicalLaw,
        nonLinearCowperSymondsPlastic,
        nonLinGeomMechLaw,
        nonLinearCowperSymondsPlastic
    );
}

// ************************************************************************* //
