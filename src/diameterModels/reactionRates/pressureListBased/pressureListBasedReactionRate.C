/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "pressureListBasedReactionRate.H"
#include "thermodynamicConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRates
{
    defineTypeNameAndDebug(pressureListBased, 0);
    addToRunTimeSelectionTable(reactionRate, pressureListBased, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reactionRates::pressureListBased::pressureListBased(const dictionary& dict)
:
    reactionRate(dict),
    pScale_(dict.lookup<scalar>("pScale")),
    pExponent_(dict.lookup<scalarList>("pExponent")),
    pCoeff_(dict.lookup<scalarList>("pCoeff")),
    pRange_(dict.lookup<scalarList>("pRange"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reactionRates::pressureListBased::~pressureListBased()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::reactionRates::pressureListBased::k
(
    const scalar& p,
    const scalar& T
) const
{
    //Find what the pressure value is and set coefficients accordingly
    // pRange = [100 200 300]
    // p.value = 50 -> alpha(1), beta(1)
    // p.value = 150 ->alpha(2), beta(2)
    // p.value = 350 ->alpha(3), beta(3)
    label i;
   if (pRange_.size() == 1)
    {
        i = 0;
    }
    else if (p > pRange_[pRange_.size()-1])
    {
        i = pRange_.size()-1;
    } else
    {
        for(i = 1; i<pRange_.size(); ++i)
        {
            if(p < pRange_[0] || p == pRange_[i-1])
            {
                --i;
                break;
            }
            else if(p > pRange_[i-1] && p < pRange_[i])
            {
                break;
            }        
        }
    }

    return pCoeff_[i]*pow(p*pScale_, pExponent_[i]);
}


Foam::tmp<Foam::volScalarField> Foam::reactionRates::pressureListBased::k
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> tmpk
    (
        new volScalarField
        (
            IOobject
            (
                "pressureListBased:k",
                p.time().timeName(),
                p.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            p.mesh(),
            dimensionedScalar(dimVelocity,0)
        )
    );
    volScalarField& K = tmpk.ref();
    forAll (K,celli)
    {
       K[celli] = k(p[celli],T[celli]);
    }
    return tmpk;
}

// ************************************************************************* //
