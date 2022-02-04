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

#include "JohnsonCookPlastic.H"
#include "fvc.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class PlasticType>
Foam::scalar Foam::JohnsonCookPlastic<PlasticType>::curYieldStress
(
    const scalar epsilonPEqOld,     // Old equivalent plastic strain
    const scalar curEpsilonPEq,     // Current equivalent plastic strain
    const scalar J                  // Current Jacobian
) const
{
    scalar eDot
    (
        curEpsilonPEq - epsilonPEqOld/this->mesh().time().deltaTValue()
    );
    if (totalStrainRate_)
    {
        eDot += eDoti_;
    }
    scalar sigmay =
        (
            A_.value()
          + B_.value()
           *pow(max(curEpsilonPEq/epsilonP0_.value(), 0), n_.value())
        )
       *(
            1.0
          + (C_.value()*log(max(eDot/epsilonPDot0_.value(), 1.0)))
        );
    if (temperatureEffects_)
    {
        sigmay *= (1.0 - (Ti_ - Tref_.value())/(Tm_.value() - Tref_.value()));
    }
    return J*sigmay;
}


template<class PlasticType>
void Foam::JohnsonCookPlastic<PlasticType>::setCellValues(const label celli) const
{
    Ti_ = T_()[celli];
    eDoti_ = epsilonElasticEffDot_()[celli];
}


template<class PlasticType>
void Foam::JohnsonCookPlastic<PlasticType>::setVolPatchFaceValues
(
    const label patchi,
    const label facei
) const
{
    Ti_ = T_().boundaryField()[patchi][facei];
    eDoti_ = epsilonElasticEffDot_().boundaryField()[patchi][facei];
}


template<class PlasticType>
void Foam::JohnsonCookPlastic<PlasticType>::setFaceValues(const label facei) const
{
    Ti_ = Tf_()[facei];
    eDoti_ = epsilonElasticEffDotf_()[facei];
}


template<class PlasticType>
void Foam::JohnsonCookPlastic<PlasticType>::setSurfacePatchFaceValues
(
    const label patchi,
    const label facei
) const
{
    Ti_ = Tf_().boundaryField()[patchi][facei];
    eDoti_ = epsilonElasticEffDotf_().boundaryField()[patchi][facei];
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class PlasticType>
Foam::JohnsonCookPlastic<PlasticType>::JohnsonCookPlastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    PlasticType(name, mesh, dict, nonLinGeom),
    T_(),
    A_("A", dimPressure, dict),
    B_("B", dimPressure, dict),
    C_("C", dimless, dict),
    n_("n", dimless, dict),
    epsilonP0_
    (
        dimensionedScalar::lookupOrDefault
        (
            "epsilonP0",
            dict,
            dimless,
            1.0
        )
    ),
    epsilonPDot0_
    (
        dimensionedScalar::lookupOrDefault
        (
            "epsilonPDot0",
            dict,
            inv(dimTime),
            1.0
        )
    ),
    totalStrainRate_(dict.lookupOrDefault<Switch>("totalStrainRate", false)),
    temperatureEffects_
    (
        dict.lookupOrDefault<Switch>("temperatureEffects", false)
    ),
    m_(1.0),
    Tref_("Tref", dimTemperature, 300),
    Tm_("Tm", dimTemperature, 1000)
{
    if (temperatureEffects_)
    {
        m_ = dict.lookup<scalar>("m");
        Tref_.read(dict);
        Tm_.read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class PlasticType>
Foam::JohnsonCookPlastic<PlasticType>::~JohnsonCookPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class PlasticType>
void Foam::JohnsonCookPlastic<PlasticType>::correct
(
    volSymmTensorField& sigma
)
{
    if (!T_.valid())
    {
        T_.set(&this->mesh().template lookupObject<volScalarField>("T"));
    }

    // Compute effective strain rate
    this->updateEpsilon
    (
        this->epsilonRef(),
        this->nu_/this->E_,
        sigma,
        this->epsilonP()
    );
    tmp<volSymmTensorField> eDot(fvc::ddt(this->epsilon()));
    const volSymmTensorField eElasticDot
    (
        dev(eDot) - dev(fvc::ddt(this->epsilonP()))
    );
    epsilonElasticEffDot_ = sqrt(2.0/3.0*(eElasticDot && eElasticDot));

    PlasticType::correct(sigma, false);

    epsilonElasticEffDot_.clear();
}


template<class PlasticType>
void Foam::JohnsonCookPlastic<PlasticType>::correct
(
    surfaceSymmTensorField& sigma
)
{
    if (!T_.valid())
    {
        T_.set(&this->mesh().template lookupObject<volScalarField>("T"));
    }
    Tf_ = fvc::interpolate(T_());

    // Compute effective strain rate
    this->updateEpsilon
    (
        this->epsilonfRef(),
        this->nu_/this->E_,
        sigma,
        this->epsilonPf()
    );
    tmp<surfaceSymmTensorField> eDot(symm(fvc::ddt(this->epsilonf())));
    const surfaceSymmTensorField eElasticDot
    (
        dev(eDot) - dev(fvc::ddt(this->epsilonPf()))
    );
    epsilonElasticEffDotf_ = sqrt(2.0/3.0*(eElasticDot && eElasticDot));

    PlasticType::correct(sigma, false);

    epsilonElasticEffDotf_.clear();
}

// ************************************************************************* //
