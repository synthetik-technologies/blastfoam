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

#include "CowperSymondsPlastic.H"
#include "fvc.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class PlasticType>
Foam::scalar Foam::CowperSymondsPlastic<PlasticType>::curYieldStress
(
    const scalar epsilonPEqOld,     // Old equivalent plastic strain
    const scalar curEpsilonPEq,     // Current equivalent plastic strain
    const scalar J                  // Current Jacobian
) const
{
    scalar eDot
    (
        mag(curEpsilonPEq - epsilonPEqOld/this->mesh().time().deltaTValue())
    );
    if (totalStrainRate_)
    {
        eDot += eDoti_;
    }
    scalar sigmay =
//         (
//             this->K_.value()
//            *pow
//             (
//                 curEpsilonPEq
//               + pow(this->E_.value()/this->K_.value(), 1/(n_ - 1.0)),
//                 n_
//             )
//         )
        sigmaYi_
       *pow(1.0 + eDot/C_.value(), 1.0/p_);
    if (temperatureEffects_)
    {
        sigmay *= (1.0 - (Ti_ - Tref_.value())/(Tm_.value() - Tref_.value()));
    }
    return sigmay;
}


template<class PlasticType>
void Foam::CowperSymondsPlastic<PlasticType>::setCellValues(const label celli) const
{
    Ti_ = T_()[celli];
    eDoti_ = epsilonElasticEffDotPtr_()[celli];
    sigmaYi_ = this->sigmaYPtr_()[celli];
}


template<class PlasticType>
void Foam::CowperSymondsPlastic<PlasticType>::setVolPatchFaceValues
(
    const label patchi,
    const label facei
) const
{
    Ti_ = T_().boundaryField()[patchi][facei];
    eDoti_ = epsilonElasticEffDotPtr_().boundaryField()[patchi][facei];
    sigmaYi_ = this->sigmaYPtr_().boundaryField()[patchi][facei];
}


template<class PlasticType>
void Foam::CowperSymondsPlastic<PlasticType>::setFaceValues(const label facei) const
{
    Ti_ = Tf_()[facei];
    eDoti_ = epsilonElasticEffDotfPtr_()[facei];
    sigmaYi_ = this->sigmaYfPtr_()[facei];
}


template<class PlasticType>
void Foam::CowperSymondsPlastic<PlasticType>::setSurfacePatchFaceValues
(
    const label patchi,
    const label facei
) const
{
    Ti_ = Tf_().boundaryField()[patchi][facei];
    eDoti_ = epsilonElasticEffDotfPtr_().boundaryField()[patchi][facei];
    sigmaYi_ = this->sigmaYfPtr_().boundaryField()[patchi][facei];
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class PlasticType>
Foam::CowperSymondsPlastic<PlasticType>::CowperSymondsPlastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    PlasticType(name, mesh, dict, nonLinGeom),
    T_(),
    C_("C", dimless, dict),
    p_(dict.lookup<scalar>("p")),
    n_(dict.lookup<scalar>("n")),
    sigmaY0_
    (
        "sigmaY",
        this->K_*pow(pow(this->E_/this->K_, 1.0/(n_ - 1.0)), n_)
    ),
    totalStrainRate_(dict.lookupOrDefault<Switch>("totalStrainRate", true)),
    temperatureEffects_
    (
        dict.lookupOrDefault<Switch>("temperatureEffects", true)
    ),
    m_(1.0),
    Tref_("Tref", dimTemperature, 300),
    Tm_("Tm", dimTemperature, 1000)
{
    sigmaY0_.readIfPresent(dict);
    if (temperatureEffects_)
    {
        m_ = dict.lookup<scalar>("m");
        Tref_.read(dict);
        Tm_.read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class PlasticType>
Foam::CowperSymondsPlastic<PlasticType>::~CowperSymondsPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class PlasticType>
void Foam::CowperSymondsPlastic<PlasticType>::correct
(
    volSymmTensorField& sigma
)
{
    if (!T_.valid())
    {
        T_.set(&this->mesh().template lookupObject<volScalarField>("T"));
    }

    // Compute effective strain rate
    tmp<volSymmTensorField> eDot(fvc::ddt(this->epsilon()));
    const volSymmTensorField eElasticDot
    (
        dev(eDot) - dev(fvc::ddt(this->epsilonP()))
    );
    epsilonElasticEffDotRef() = sqrt(2.0/3.0*(eElasticDot && eElasticDot));

    PlasticType::correct(sigma);

}


template<class PlasticType>
void Foam::CowperSymondsPlastic<PlasticType>::correct
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
    const surfaceTensorField& gradD =
        this->mesh().template lookupObject<surfaceTensorField>("grad(D)f");
    tmp<surfaceSymmTensorField> eDot(symm(fvc::ddt(gradD)));
    const surfaceSymmTensorField eElasticDot
    (
        dev(eDot) - dev(fvc::ddt(this->epsilonPf()))
    );
    epsilonElasticEffDotfRef() = sqrt(2.0/3.0*(eElasticDot && eElasticDot));

    PlasticType::correct(sigma);
}

// ************************************************************************* //
