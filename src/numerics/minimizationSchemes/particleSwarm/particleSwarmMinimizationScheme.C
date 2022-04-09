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

#include "particleSwarmMinimizationScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(particleSwarmMinimizationScheme, 0);
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        particleSwarmMinimizationScheme,
        dictionaryUnivariate
    );
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        particleSwarmMinimizationScheme,
        dictionaryMultivariate
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleSwarmMinimizationScheme::particleSwarmMinimizationScheme
(
    const scalarUnivariateEquation& eqns,
    const dictionary& dict
)
:
    minimizationScheme(eqns, dict),
    rand_(0),
    particles_(dict.lookupOrDefault<label>("nParticles", 100)),
    cLocal_(dict.lookup<scalar>("cLocal")),
    cGlobal_(dict.lookup<scalar>("cGlobal")),
    vWeight_(dict.lookup<scalar>("vWeight"))

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::particleSwarmMinimizationScheme::minimize
(
    const scalarField& x0,
    const scalarField& xLow,
    const scalarField& xHigh,
    const label li
) const
{
    label n = eqns_.nVar();
    label np = particles_.size();
    tmp<scalarField> txMean(new scalarField(n, 0.0));
    scalarField& xMean = txMean.ref();

    scalar yBest = great;
    scalar y = yBest;
    scalarField xBest(xMean);
    scalarField xVar(n, 0.0);
    forAll(particles_, i)
    {
        particle& p = particles_[i];
        p.x.resize(n);
        for (label j = 0; j < n; j++)
        {
            p.x[j] =
                rand_.sampleAB<scalar>(xLow[j], xHigh[j]);
        }
        p.xBest = p.x;
        p.v.resize(n, 0.0);

        y = eqns_.fX(p.x, li);
        p.y = y;
        p.yBest = p.y;
        if (y < yBest)
        {
            yBest = y;
            xBest = p.x;
        }
        xMean += p.x;
    }
    xMean /= scalar(np);
    forAll(particles_, i)
    {
        xVar += sqr(particles_[i].x - xMean);
    }
    xVar /= scalar(np);
    scalarField xStd(sqrt(xVar));
    if (normalize_)
    {
        xStd /= stabilise(xMean, small);
    }

    scalarField r1(n);
    scalarField r2(n);
    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        if (converged(xStd))
        {
            break;
        }

        xMean = 0;
        forAll(particles_, i)
        {
            particle& p = particles_[i];
            forAll(r1, cmpti)
            {
                r1[cmpti] = rand_.sample01<scalar>();
                r2[cmpti] = rand_.sample01<scalar>();
            }
            scalarField xOld(p.x);
            p.x += p.v;
            eqns_.limit(p.x);
            p.v = p.x - xOld;
            p.v =
                vWeight_*p.v
              + cLocal_*r1*(p.xBest - p.x)
              + cGlobal_*r2*(xBest - p.x);
            y = eqns_.fX(p.x, li);
            if (y < yBest)
            {
                yBest = y;
                xBest = p.xBest;
            }
            if (y < p.yBest)
            {
                p.yBest = y;
                p.xBest = p.x;
            }
            xMean += p.x;
        }
        xMean /= scalar(np);

        xVar = 0.0;
        forAll(particles_, i)
        {
            xVar += sqr(particles_[i].x - xMean);
        }
        xVar /= scalar(np);
        xStd = sqrt(xVar);
        if (normalize_)
        {
            xStd /= stabilise(xMean, small);
        }

        printStepInformation(xMean);
    }
    xMean = xBest;
    printFinalInformation();
    return txMean;
}

// ************************************************************************* //
