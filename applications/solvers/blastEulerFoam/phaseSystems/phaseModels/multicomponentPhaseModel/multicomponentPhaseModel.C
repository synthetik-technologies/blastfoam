/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2019-04-29 Jeff Heylmun:    Simplified model
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "multicomponentPhaseModel.H"
#include "phaseSystem.H"
#include "fvMatrix.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multicomponentPhaseModel, 0);
    addToRunTimeSelectionTable
    (
        phaseModel,
        multicomponentPhaseModel,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multicomponentPhaseModel::multicomponentPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    phaseModel(fluid, phaseName, index),
    components_(phaseDict_.lookup("components")),
    p_
    (
        IOobject
        (
            IOobject::groupName("p", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid_.p(),
        fluid_.p().boundaryField()
    ),
    rho_
    (
        IOobject
        (
            IOobject::groupName("rho", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", name_),
            this->mesh().time().timeName(),
            this->mesh()
        ),
        this->mesh(),
        dimensionedScalar("0", phi_.dimensions(), 0.0)
    ),
    thermos_(components_.size()),
    Ys_(components_.size()),
    residualAlpha_("residualAlpha", dimless, phaseDict_),
    alphaRhoYs_(components_.size()),
    alphaRhoYPhis_(components_.size()),
    alphaRhoYsOld_(components_.size()),
    deltaAlphaRhoYs_(components_.size()),
    fluxScheme_(fluxScheme::New(fluid.mesh(), name_))
{
    //- Temporarily Store read density
    volScalarField rhoOld(rho_);
    volScalarField sumY
    (
        IOobject
        (
            "sumY",
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        0.0
    );

    forAll(components_, i)
    {
        Ys_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Y", components_[i]),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                rho_.mesh()
            )
        );
        sumY += Ys_[i];

        alphaRhoYs_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaRho", components_[i]),
                    fluid.mesh().time().timeName(),
                    fluid.mesh()
                ),
                (*this)*Ys_[i]*rho_
            )
        );
        alphaRhoYPhis_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaRhoPhi", components_[i]),
                    fluid.mesh().time().timeName(),
                    fluid.mesh()
                ),
                fluid.mesh(),
                dimensionedScalar(dimDensity*dimVelocity*dimArea, 0.0)
            )
        );

        thermos_.set
        (
            i,
            fluidThermoModel::New
            (
                components_[i],
                p_,
                rho_,
                e_,
                T_,
                phaseDict_.subDict(components_[i]),
                false,
                phaseName
            ).ptr()
        );

        alphaRhoYsOld_.set(i, new PtrList<volScalarField>());
        deltaAlphaRhoYs_.set(i, new PtrList<volScalarField>());
    }
    // Reset density to correct value
    rho_ = rhoOld;

    if (Foam::min(sumY).value() != 1 && Foam::max(sumY).value() != 1)
    {
        FatalErrorInFunction
            << "The sum of components is not equal to 1" << nl
            << abort(FatalError);
    }

    volScalarField e("e", Ys_[0]*thermos_[0].calce());
    for (label i = 1; i < Ys_.size(); i++)
    {
        e += Ys_[i]*thermos_[i].calce();
    }
    e_ = e;
    forAll(e_.boundaryField(), patchi)
    {
        forAll(e_.boundaryField()[patchi], facei)
        {
            e_.boundaryFieldRef()[patchi][facei] =
                e.boundaryField()[patchi][facei];
        }
    }
    e_.correctBoundaryConditions();


    T_ = Ys_[0]*thermos_[0].calcT();
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        T_ += Ys_[phasei]*thermos_[phasei].calcT();
    }
    T_.correctBoundaryConditions();
    correctThermo();

    if (Foam::max(mu()).value() > 0)
    {
        turbulence_ =
            PhaseCompressibleTurbulenceModel<phaseModel>::New
            (
                *this,
                rho_,
                U_,
                alphaRhoPhi_,
                phi_,
                *this
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multicomponentPhaseModel::~multicomponentPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multicomponentPhaseModel::solveAlpha(const bool s)
{
    phaseModel::solveAlpha(s);
}


void Foam::multicomponentPhaseModel::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    PtrList<volScalarField> alphaRhoYsOld(Ys_.size());
    forAll(Ys_, i)
    {
        alphaRhoYsOld.set
        (
            i, new volScalarField(alphaRhoYs_[i])
        );
        this->storeOld(stepi, alphaRhoYsOld[i], alphaRhoYsOld_[i]);
        this->blendOld(stepi, alphaRhoYsOld[i], alphaRhoYsOld_[i], ai);
    }

    PtrList<volScalarField> deltaAlphaRhoYs(Ys_.size());
    forAll(Ys_, i)
    {
        deltaAlphaRhoYs.set
        (
            i,
            new volScalarField(fvc::div(alphaRhoYPhis_[i]))
        );
        this->storeDelta
        (
            stepi,
            deltaAlphaRhoYs[i],
            deltaAlphaRhoYs_[i]
        );
        this->blendDelta
        (
            stepi,
            deltaAlphaRhoYs[i],
            deltaAlphaRhoYs_[i],
            bi
        );
    }

    dimensionedScalar dT = rho_.time().deltaT();

    forAll(Ys_, i)
    {
        alphaRhoYs_[i].oldTime() = alphaRhoYsOld[i];
        alphaRhoYs_[i] =
            alphaRhoYsOld[i] - dT*deltaAlphaRhoYs[i];
        alphaRhoYs_[i].correctBoundaryConditions();

        thermos_[i].solve(stepi, ai, bi);
    }

    if (solveAlpha_)
    {
        volScalarField& alpha = *this;
        volScalarField alphaOld(alpha);
        this->storeOld(stepi, alphaOld, alphaOld_);
        this->blendOld(stepi, alphaOld, alphaOld_, ai);

        volScalarField deltaAlpha
        (
            (fluid_.U() & fvc::grad(fluxScheme_->alphaf()))
        );

        this->storeDelta(stepi, deltaAlpha, deltaAlpha_);
        this->blendDelta(stepi, deltaAlpha, deltaAlpha_, bi);


        dimensionedScalar dT = rho_.time().deltaT();

        alpha = alphaOld - dT*(deltaAlpha);
        alpha.max(0);
        alpha.min(alphaMax_);
        alpha.correctBoundaryConditions();
    }

    phaseModel::solve(stepi, ai, bi);
}


void Foam::multicomponentPhaseModel::postUpdate()
{
    phaseModel::postUpdate();

    forAll(thermos_, i)
    {
        thermos_[i].postUpdate();
    }
}


void Foam::multicomponentPhaseModel::clearODEFields()
{
    phaseModel::clearODEFields();
    fluxScheme_->clear();

    if (solveAlpha_)
    {
        this->clearOld(alphaOld_);
        this->clearDelta(deltaAlpha_);
    }

    forAll(alphaRhoYsOld_, i)
    {
        this->clearOld(alphaRhoYsOld_[i]);
        this->clearDelta(deltaAlphaRhoYs_[i]);
        thermos_[i].clearODEFields();
    }
}


void Foam::multicomponentPhaseModel::update()
{
    tmp<volScalarField> c(speedOfSound());
    const volScalarField& alpha = *this;

    fluxScheme_->update
    (
        alpha,
        rho_,
        U_,
        e_,
        p_,
        c(),
        phi_,
        alphaPhi_,
        alphaRhoPhi_,
        alphaRhoUPhi_,
        alphaRhoEPhi_
    );
    forAll(Ys_, i)
    {
        alphaRhoYPhis_[i] =
            alphaRhoPhi_*fluxScheme_->interpolate(Ys_[i], Ys_[i].name());
    }
}


void Foam::multicomponentPhaseModel::decode()
{
    this->correctBoundaryConditions();
    volScalarField alpha(Foam::max(*this, residualAlpha()));

    alphaRho_.max(0);
    rho_.ref() = alphaRho_/alpha();
    rho_.correctBoundaryConditions();
    alphaRho_.boundaryFieldRef() = this->boundaryField()*rho_.boundaryField();

    volScalarField alphaRho(alphaRho_);
    alphaRho.max(1e-10);
    volScalarField sumY
    (
        IOobject
        (
            "sumY",
            rho_.time().timeName(),
            rho_.mesh()
        ),
        rho_.mesh(),
        0.0
    );
    volScalarField sumRhoY
    (
        IOobject
        (
            "sumRhoY",
            rho_.time().timeName(),
            rho_.mesh()
        ),
        rho_.mesh(),
        dimensionedScalar(dimDensity, 0.0)
    );
    for (label i = 0; i < Ys_.size() - 1; i++)
    {
        alphaRhoYs_[i].max(0);
        Ys_[i] = alphaRhoYs_[i]/alphaRho;
        Ys_[i].min(1);
        Ys_[i].correctBoundaryConditions();
        sumY += Ys_[i];
    }

    label i = Ys_.size() - 1;
    Ys_[i] = 1.0 - sumY;
    Ys_[i].max(0.0);
    Ys_[i].min(1.0);
    Ys_[i].correctBoundaryConditions();
    alphaRhoYs_[i] = (*this)*Ys_[i]*rho_;

    U_.ref() = alphaRhoU_()/(alphaRho());
    U_.correctBoundaryConditions();

    alphaRhoU_.boundaryFieldRef() =
        (*this).boundaryField()*rho_.boundaryField()*U_.boundaryField();

    volScalarField E(alphaRhoE_/alphaRho);
    e_.ref() = E() - 0.5*magSqr(U_());

//     //--- Hard limit, e
//     if(Foam::min(e_).value() < 0)
//     {
//         e_.max(small);
//         alphaRhoE_.ref() = (*this)()*rho_*(e_() + 0.5*magSqr(U_()));
//     }
    e_.correctBoundaryConditions();

    alphaRhoE_.boundaryFieldRef() =
        (*this).boundaryField()
       *rho_.boundaryField()
       *(
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );

    correctThermo();
}


void Foam::multicomponentPhaseModel::encode()
{
    for (label i = 0; i < Ys_.size(); i++)
    {
        alphaRhoYs_[i] = (*this)*rho_*Ys_[i];
    }

    alphaRho_ = (*this)*rho_;
    alphaRhoU_ = alphaRho_*U_;
    alphaRhoE_ = alphaRho_*(e_ + 0.5*magSqr(U_));
}


Foam::tmp<Foam::volVectorField>
Foam::multicomponentPhaseModel::gradP() const
{
    return fvc::grad(fluxScheme_->pf());
}


Foam::tmp<Foam::volVectorField>
Foam::multicomponentPhaseModel::gradAlpha() const
{
    return fvc::grad(fluxScheme_->alphaf());
}



void Foam::multicomponentPhaseModel::correctThermo()
{
    volScalarField pByGamma
    (
        IOobject
        (
            "pRef",
            e_.mesh().time().timeName(),
            e_.mesh()
        ),
        e_.mesh(),
        dimensionedScalar("0", dimPressure, 0.0)
    );
    volScalarField rGamma
    (
        IOobject
        (
            "rGamma",
            e_.mesh().time().timeName(),
            e_.mesh()
        ),
        e_.mesh(),
        0.0
    );

    forAll(Ys_, phasei)
    {
        volScalarField alphaByGamma
        (
            Ys_[phasei]/(thermos_[phasei].Gamma() - 1.0)
        );
        rGamma += alphaByGamma;
        pByGamma += thermos_[phasei].calcP()*alphaByGamma;
    }
    rGamma.max(0.0);

    p_ = pByGamma/rGamma;
    p_.max(small);
    p_.correctBoundaryConditions();

    T_ = Ys_[0]*thermos_[0].calcT();
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        T_ += Ys_[phasei]*thermos_[phasei].calcT();
    }
    T_.correctBoundaryConditions();

    forAll(thermos_, i)
    {
        thermos_[i].correct();
    }
}


Foam::tmp<Foam::volScalarField>
Foam::multicomponentPhaseModel::speedOfSound() const
{
    volScalarField alphaXiRhoCSqr
    (
        IOobject
        (
            "alphaXiRhoCSqr",
            e_.mesh().time().timeName(),
            e_.mesh()
        ),
        e_.mesh(),
        dimensionedScalar("0", dimPressure, 0.0)
    );
    volScalarField Xi
    (
        IOobject
        (
            "Xi",
            e_.mesh().time().timeName(),
            e_.mesh()
        ),
        e_.mesh(),
        0.0
    );

    volScalarField rho(rho_);
    rho.max(1e-10);
    forAll(thermos_, phasei)
    {
        tmp<volScalarField> alphaXi
        (
            Ys_[phasei]/(thermos_[phasei].Gamma() - 1.0)
        );

        Xi += alphaXi();
        alphaXiRhoCSqr += alphaXi*rho_*sqr(thermos_[phasei].speedOfSound());
    }

    tmp<volScalarField> cSqr(alphaXiRhoCSqr/(rho*Xi));
    cSqr.ref().max(small);
    return sqrt(cSqr);
}


Foam::tmp<Foam::volScalarField>
Foam::multicomponentPhaseModel::ESource() const
{
    tmp<volScalarField> ESourceTmp(Ys_[0]*thermos_[0].ESource());
    for (label i = 1; i < thermos_.size(); i++)
    {
        ESourceTmp.ref() += Ys_[i]*thermos_[i].ESource();
    }
    return ESourceTmp*(*this);
}


Foam::tmp<Foam::volScalarField> Foam::multicomponentPhaseModel::Cv() const
{
    tmp<volScalarField> CvTmp(Ys_[0]*thermos_[0].Cv());
    for (label i = 1; i < thermos_.size(); i++)
    {
        CvTmp.ref() += Ys_[i]*thermos_[i].Cv();
    }
    return CvTmp;
}

Foam::tmp<Foam::volScalarField> Foam::multicomponentPhaseModel::Cp() const
{
    tmp<volScalarField> CpTmp(Ys_[0]*thermos_[0].Cp());
    for (label i = 1; i < thermos_.size(); i++)
    {
        CpTmp.ref() += Ys_[i]*thermos_[i].Cp();
    }
    return CpTmp;
}

Foam::scalar Foam::multicomponentPhaseModel::Cvi(const label celli) const
{
    scalar cv(Ys_[0][celli]*thermos_[0].Cvi(celli));
    for (label i = 1; i < thermos_.size(); i++)
    {
        cv += Ys_[i][celli]*thermos_[i].Cvi(celli);
    }
    return cv;
}

Foam::tmp<Foam::volScalarField> Foam::multicomponentPhaseModel::mu() const
{
    tmp<volScalarField> muTmp(Ys_[0]*thermos_[0].mu());
    for (label i = 1; i < thermos_.size(); i++)
    {
        muTmp.ref() += Ys_[i]*thermos_[i].mu();
    }
    return muTmp;
}


Foam::tmp<Foam::scalarField>
Foam::multicomponentPhaseModel::mu(const label patchi) const
{
    tmp<scalarField> muTmp
    (
        Ys_[0].boundaryField()[patchi]*thermos_[0].mu(patchi)
    );
    for (label i = 1; i < thermos_.size(); i++)
    {
        muTmp.ref() +=
            Ys_[i].boundaryField()[patchi]*thermos_[i].mu(patchi);
    }
    return muTmp;
}

Foam::tmp<Foam::volScalarField> Foam::multicomponentPhaseModel::nu() const
{
    tmp<volScalarField> nuTmp(Ys_[0]*thermos_[0].nu());
    for (label i = 1; i < thermos_.size(); i++)
    {
        nuTmp.ref() += Ys_[i]*thermos_[i].nu();
    }
    return nuTmp;
}

Foam::tmp<Foam::scalarField>
Foam::multicomponentPhaseModel::nu(const label patchi) const
{
    tmp<scalarField> nuTmp
    (
        Ys_[0].boundaryField()[patchi]*thermos_[0].nu(patchi)
    );
    for (label i = 1; i < thermos_.size(); i++)
    {
        nuTmp.ref() +=
            Ys_[i].boundaryField()[patchi]*thermos_[i].nu(patchi);
    }
    return nuTmp;
}

Foam::scalar
Foam::multicomponentPhaseModel::nui(const label celli) const
{
    scalar n(Ys_[0][celli]*thermos_[0].nui(celli));
    for (label i = 1; i < thermos_.size(); i++)
    {
        n += Ys_[i][celli]*thermos_[i].nui(celli);
    }
    return n;
}

Foam::tmp<Foam::volScalarField> Foam::multicomponentPhaseModel::alpha() const
{
    tmp<volScalarField> alphaTmp(Ys_[0]*thermos_[0].alpha());
    for (label i = 1; i < thermos_.size(); i++)
    {
        alphaTmp.ref() += Ys_[i]*thermos_[i].alpha();
    }
    return alphaTmp;
}

Foam::tmp<Foam::scalarField>
Foam::multicomponentPhaseModel::alpha(const label patchi) const
{
    tmp<scalarField> alphaTmp
    (
        Ys_[0].boundaryField()[patchi]*thermos_[0].alpha(patchi)
    );
    for (label i = 1; i < thermos_.size(); i++)
    {
        alphaTmp.ref() +=
            Ys_[i].boundaryField()[patchi]*thermos_[i].alpha(patchi);
    }
    return alphaTmp;
}

Foam::tmp<Foam::volScalarField> Foam::multicomponentPhaseModel::alphaEff
(
    const volScalarField& alphat
) const
{
    tmp<volScalarField> alphaEffTmp(Ys_[0]*thermos_[0].alphaEff(alphat));
    for (label i = 1; i < thermos_.size(); i++)
    {
        alphaEffTmp.ref() += Ys_[i]*thermos_[i].alphaEff(alphat);
    }
    return alphaEffTmp;
}

Foam::tmp<Foam::scalarField> Foam::multicomponentPhaseModel::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    tmp<scalarField> alphaEffTmp
    (
        Ys_[0].boundaryField()[patchi]
       *thermos_[0].alphaEff(alphat, patchi)
    );
    for (label i = 1; i < thermos_.size(); i++)
    {
        alphaEffTmp.ref() +=
            Ys_[i].boundaryField()[patchi]
           *thermos_[i].alphaEff(alphat, patchi);
    }
    return alphaEffTmp;
}

Foam::tmp<Foam::volScalarField> Foam::multicomponentPhaseModel::alphahe() const
{
    tmp<volScalarField> alphaheTmp(Ys_[0]*thermos_[0].alphahe());
    for (label i = 1; i < thermos_.size(); i++)
    {
        alphaheTmp.ref() += Ys_[i]*thermos_[i].alphahe();
    }
    return alphaheTmp;
}

Foam::tmp<Foam::scalarField>
Foam::multicomponentPhaseModel::alphahe(const label patchi) const
{
    tmp<scalarField> alphaheTmp
    (
        Ys_[0].boundaryField()[patchi]*thermos_[0].alphahe(patchi)
    );
    for (label i = 1; i < thermos_.size(); i++)
    {
        alphaheTmp.ref() +=
            Ys_[i].boundaryField()[patchi]*thermos_[i].alphahe(patchi);
    }
    return alphaheTmp;
}

Foam::tmp<Foam::volScalarField> Foam::multicomponentPhaseModel::kappa() const
{
    tmp<volScalarField> kappaTmp(Ys_[0]*thermos_[0].kappa());
    for (label i = 1; i < thermos_.size(); i++)
    {
        kappaTmp.ref() += Ys_[i]*thermos_[i].kappa();
    }
    return kappaTmp;
}

Foam::tmp<Foam::scalarField>
Foam::multicomponentPhaseModel::kappa(const label patchi) const
{
    tmp<scalarField> kappaTmp
    (
        Ys_[0].boundaryField()[patchi]*thermos_[0].kappa(patchi)
    );
    for (label i = 1; i < thermos_.size(); i++)
    {
        kappaTmp.ref() +=
            Ys_[i].boundaryField()[patchi]*thermos_[i].kappa(patchi);
    }
    return kappaTmp;
}

Foam::scalar
Foam::multicomponentPhaseModel::kappai(const label celli) const
{
    scalar k(Ys_[0][celli]*thermos_[0].kappai(celli));
    for (label i = 1; i < thermos_.size(); i++)
    {
        k += Ys_[i][celli]*thermos_[i].kappai(celli);
    }
    return k;
}

Foam::tmp<Foam::volScalarField> Foam::multicomponentPhaseModel::kappaEff
(
    const volScalarField& alphat
) const
{
    tmp<volScalarField> kappaEffTmp(Ys_[0]*thermos_[0].kappaEff(alphat));
    for (label i = 1; i < thermos_.size(); i++)
    {
        kappaEffTmp.ref() +=
            Ys_[i]*thermos_[i].kappaEff(alphat);
    }
    return kappaEffTmp;
}

Foam::tmp<Foam::scalarField>
Foam::multicomponentPhaseModel::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    tmp<scalarField> kappaEffTmp
    (
        Ys_[0].boundaryField()[patchi]*thermos_[0].kappaEff(alphat, patchi)
    );
    for (label i = 1; i < thermos_.size(); i++)
    {
        kappaEffTmp.ref() +=
            Ys_[i].boundaryField()[patchi]*thermos_[i].kappaEff(alphat, patchi);
    }
    return kappaEffTmp;
}

// ************************************************************************* //
