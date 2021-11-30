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
          + C_.value()*log(max(eDot/epsilonPDot0_.value(), 1e-6))
        );
    if (temperatureEffects_)
    {
        sigmay *= (1.0 - (Ti_ - Tref_.value())/(Tm_.value() - Tref_.value()));
    }
    return sigmay;
}


template<class PlasticType>
void Foam::JohnsonCookPlastic<PlasticType>::setCellValues(const label celli) const
{
    Ti_ = T_[celli];
    eDoti_ = epsilonElasticEffDot_()[celli];
}


template<class PlasticType>
void Foam::JohnsonCookPlastic<PlasticType>::setVolPatchFaceValues
(
    const label patchi,
    const label facei
) const
{
    Ti_ = T_.boundaryField()[patchi][facei];
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
    T_(mesh.lookupObject<volScalarField>("T")),
    A_("A", dimPressure, dict),
    B_("B", dimPressure, dict),
    C_("C", dimless, dict),
    n_("n", dimless, dict),
    m_("m", dimless, dict),
    Tref_("Tref", dimTemperature, dict),
    Tm_("Tm", dimTemperature, dict),
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
    totalStrainRate_(dict.lookupOrDefault<Switch>("totalStrainRate", true)),
    temperatureEffects_
    (
        dict.lookupOrDefault<Switch>("temperatureEffects", true)
    )
{
    this->sigmaY_ = A_;
    this->sigmaYf_ = A_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class PlasticType>
Foam::JohnsonCookPlastic<PlasticType>::~JohnsonCookPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class PlasticType>
void Foam::JohnsonCookPlastic<PlasticType>::correct(volSymmTensorField& sigma)
{
    Info<<"here"<<endl;
    // Compute effective strain rate
    const volTensorField& gradD
    (
        this->mesh().template lookupObject<volTensorField>("grad(D)")
    );
    tmp<volSymmTensorField> eDot(symm(fvc::ddt(gradD)));
    const volSymmTensorField eElasticDot
    (
        dev(eDot) - dev(fvc::ddt(this->epsilonP_))
    );
    epsilonElasticEffDot_ = sqrt(2.0/3.0*(eElasticDot && eElasticDot));

    PlasticType::correct(sigma);
}


template<class PlasticType>
void Foam::JohnsonCookPlastic<PlasticType>::correct(surfaceSymmTensorField& sigma)
{
    Tf_ = fvc::interpolate(T_);

    // Compute effective strain rate
    const surfaceTensorField& gradD =
        this->mesh().template lookupObject<surfaceTensorField>("grad(D)f");
    tmp<surfaceSymmTensorField> eDot(symm(fvc::ddt(gradD)));
    const surfaceSymmTensorField eElasticDot
    (
        dev(eDot) - dev(fvc::ddt(this->epsilonPf_))
    );
    epsilonElasticEffDotf_ = sqrt(2.0/3.0*(eElasticDot && eElasticDot));

    PlasticType::correct(sigma);
}

// ************************************************************************* //
