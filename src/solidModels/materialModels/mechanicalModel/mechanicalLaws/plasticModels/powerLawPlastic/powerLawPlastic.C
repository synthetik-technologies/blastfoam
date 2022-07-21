/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
31-03-2022 Synthetik Applied Technologies: Added powerLawPlastic
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

#include "powerLawPlastic.H"
#include "fvc.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class PlasticType>
Foam::scalar Foam::powerLawPlastic<PlasticType>::curYieldStress
(
    const scalar epsilonPEqOld,     // Old equivalent plastic strain
    const scalar curEpsilonPEq,     // Current equivalent plastic strain
    const scalar J                  // Current Jacobian
) const
{
    return k_.value()*pow(curEpsilonPEq + epsilonY_.value(), n_.value());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class PlasticType>
Foam::powerLawPlastic<PlasticType>::powerLawPlastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    PlasticType(name, mesh, dict, nonLinGeom),
    k_("k", dimPressure, -1.0),
    n_("n", dimless, -1.0),
    epsilonY_("epsilonY", dimless, 0.0)
{
    scalar sigmaY(dict.lookupOrDefault<scalar>("sigmaY", -1.0));
    scalar sigmaUT(dict.lookupOrDefault<scalar>("sigmaUT", -1.0));
    scalar elongation(dict.lookupOrDefault<scalar>("elongation", -1.0));

    k_.readIfPresent(dict);
    n_.readIfPresent(dict);

    if (n_.value() > 0)
    {}
    else if (sigmaY > 0 && sigmaUT > 0 && elongation > 0)
    {
        n_.value() = log(sigmaUT/sigmaY)/log(elongation/0.002);
    }
    else
    {
        FatalErrorInFunction
            << "Strain hardening coefficient could not be determined." << nl
            << "Either provide n or sigmaY, sigmaUT, and elongation." << nl
            << abort(FatalError);
    }

    if (k_.value() > 0)
    {}
    else if (sigmaUT > 0 && sigmaY > 0)
    {
        k_.value() = sigmaY/pow(sigmaY/this->E_.value() + 0.002, n_.value());
    }
    else
    {
        FatalErrorInFunction
            << "Strength coefficient could not be determined." << nl
            << "Either provide k or sigmaY and sigmaUT." << nl
            << abort(FatalError);
    }

    if (sigmaY > 0)
    {
        epsilonY_.value() = pow(sigmaY/k_.value(), 1.0/n_.value());
    }
    else
    {
        epsilonY_.value() =
            pow(this->E_.value()/k_.value(), 1.0/(n_.value() - 1.0));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class PlasticType>
Foam::powerLawPlastic<PlasticType>::~powerLawPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
