/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

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
    k_("k", dimPressure, dict),
    n_("n", dimless, dict),
    epsilonY_("epsilonY", dimless, 0.0)
{
    if (dict.found("sigmaY"))
    {
        dimensionedScalar sigmaY("sigmaY", dimPressure, dict);
        epsilonY_ = pow(sigmaY/k_, 1.0/n_);
    }
    else
    {
        epsilonY_ = pow(this->E_/k_, 1.0/(n_ - 1.0));
    }


    this->sigmaY_ = k_*pow(epsilonY_, n_);
    this->sigmaYf_ = k_*pow(epsilonY_, n_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class PlasticType>
Foam::powerLawPlastic<PlasticType>::~powerLawPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
