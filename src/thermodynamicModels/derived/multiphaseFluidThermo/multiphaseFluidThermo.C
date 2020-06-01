/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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

#include "multiphaseFluidThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseFluidThermo::multiphaseFluidThermo
(
    const word& name,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master
)
:
    fluidThermoModel(name, p, rho, e, T, dict, master),
    phases_(dict.lookup("phases")),
    volumeFractions_(phases_.size()),
    rhos_(phases_.size()),
    thermos_(phases_.size()),
    alphaRhos_(phases_.size()),
    alphaPhis_(phases_.size()),
    alphaRhoPhis_(phases_.size())
{
    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            e.mesh().time().timeName(),
            e.mesh()
        ),
        e.mesh(),
        0.0
    );
    forAll(phases_, phasei)
    {
        word phaseName = phases_[phasei];
        volumeFractions_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alpha", phases_[phasei]),
                    rho.mesh().time().timeName(),
                    rho.mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                rho_.mesh()
            )
        );
        rhos_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("rho", phases_[phasei]),
                    rho.mesh().time().timeName(),
                    rho.mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                rho_.mesh()
            )
        );
        thermos_.set
        (
            phasei,
            fluidThermoModel::New
            (
                phaseName,
                p_,
                rhos_[phasei],
                e_,
                T_,
                dict.subDict(phaseName),
                false
            ).ptr()
        );
        thermos_[phasei].read(dict.subDict(phaseName));
        rho_ += volumeFractions_[phasei]*rhos_[phasei];
        sumAlpha += volumeFractions_[phasei];
    }
    this->initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseFluidThermo::~multiphaseFluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiphaseFluidThermo::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    forAll(phases_, phasei)
    {
        thermos_[phasei].solve(stepi, ai, bi);
    }
}


void Foam::multiphaseFluidThermo::setODEFields
(
    const label nSteps,
    const labelList& oldIs,
    const label& nOld,
    const labelList& deltaIs,
    const label nDelta
)
{
    forAll(phases_, phasei)
    {
        thermos_[phasei].setODEFields
        (
            nSteps,
            oldIs,
            nOld,
            deltaIs,
            nDelta
        );
    }
}


void Foam::multiphaseFluidThermo::clearODEFields()
{
    forAll(phases_, phasei)
    {
        thermos_[phasei].clearODEFields();
    }
}


void Foam::multiphaseFluidThermo::correct()
{
    if (master_)
    {
        T_ = calcT();
        p_ = calcP();
        p_.max(small);
    }

    forAll(thermos_, phasei)
    {
        thermos_[phasei].correct();
    }

    if (viscous_)
    {
        mu_ = volumeFractions_[0]*thermos_[0].mu();
        alpha_ = volumeFractions_[0]*thermos_[0].alpha();
        for (label phasei = 1; phasei < thermos_.size(); phasei++)
        {
            mu_ += volumeFractions_[phasei]*thermos_[phasei].mu();
            alpha_ += volumeFractions_[phasei]*thermos_[phasei].alpha();
        }
    }
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidThermo::speedOfSound() const
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

    forAll(thermos_, phasei)
    {
        tmp<volScalarField> alphaXi
        (
            volumeFractions_[phasei]/(thermos_[phasei].Gamma() - 1.0)
        );
        Xi += alphaXi();
        alphaXiRhoCSqr += alphaXi*rhos_[phasei]*sqr(thermos_[phasei].speedOfSound());
    }

    tmp<volScalarField> cSqr(alphaXiRhoCSqr/(rho_*Xi));
    cSqr.ref().max(small);
    return sqrt(cSqr);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidThermo::speedOfSound(const label patchi) const
{
    scalarField alphaXiRhoCSqr(p_.boundaryField()[patchi].size(), 0.0);
    scalarField Xi(p_.boundaryField()[patchi].size(), 0.0);

    forAll(thermos_, phasei)
    {
        scalarField alphaXi
        (
            volumeFractions_[phasei].boundaryField()[patchi]
           /(thermos_[phasei].Gamma(patchi) - 1.0)
        );
        Xi += alphaXi;
        alphaXiRhoCSqr +=
            alphaXi
           *rhos_[phasei].boundaryField()[patchi]
           *sqr(thermos_[phasei].speedOfSound(patchi));
    }

    scalarField cSqr(alphaXiRhoCSqr/(rho_.boundaryField()[patchi]*Xi));
    cSqr = max(cSqr, small);
    return sqrt(cSqr);
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseFluidThermo::calcT() const
{
    tmp<volScalarField> tmpF
    (
        volScalarField::New
        (
            IOobject::groupName("T", basicThermoModel::name_),
            volumeFractions_[0]*thermos_[0].calcT()
        )
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].calcT();
    }
    return tmpF;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidThermo::TRhoE
(
    const scalarField& T,
    const scalarField& e,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].TRhoE(T, e, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].TRhoE(T, e, patchi);
    }
    return tmpF;
}


Foam::scalar Foam::multiphaseFluidThermo::TRhoEi
(
    const scalar& T,
    const scalar& e,
    const label celli
) const
{
    scalar f(volumeFractions_[0][celli]*thermos_[0].TRhoEi(T, e, celli));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        f +=
            volumeFractions_[phasei][celli]
           *thermos_[phasei].TRhoEi(T, e, celli);
    }
    return f;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseFluidThermo::calcP() const
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

    forAll(thermos_, phasei)
    {
        volScalarField alphaByGamma
        (
            volumeFractions_[phasei]/(thermos_[phasei].Gamma() - 1.0)
        );
        rGamma += alphaByGamma;
        pByGamma += alphaByGamma*thermos_[phasei].calcP();
    }

    return pByGamma/rGamma;
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidThermo::calce() const
{
    tmp<volScalarField> tmpF
    (
        volScalarField::New
        (
            IOobject::groupName("e", basicThermoModel::name_),
            volumeFractions_[0]*thermos_[0].calce()
        )
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].calce();
    }
    return tmpF;
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidThermo::E() const
{
    tmp<volScalarField> tmpF
    (
        volScalarField::New
        (
            IOobject::groupName("E", basicThermoModel::name_),
            volumeFractions_[0]*thermos_[0].E()
        )
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].E();
    }
    return tmpF;
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidThermo::e
(
    const volScalarField& rho,
    const volScalarField& e,
    const volScalarField& T
) const
{
    tmp<volScalarField> tmpF
    (
        volScalarField::New
        (
            IOobject::groupName("E", basicThermoModel::name_),
            volumeFractions_[0]*thermos_[0].e(rho, e, T)
        )
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].e(rho, e, T);
    }
    return tmpF;
}



Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidThermo::e
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].e(rho, e, T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].e(rho, e, T, patchi);
    }
    return tmpF;
}



Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidThermo::e
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const labelList& faceCells
) const
{
    tmp<scalarField> tmpF
    (
        UIndirectList<scalar>(volumeFractions_[0], faceCells)()
       *thermos_[0].e(rho, e, T, faceCells)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            UIndirectList<scalar>(volumeFractions_[phasei], faceCells)()
           *thermos_[phasei].e(rho, e, T, faceCells);
    }
    return tmpF;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseFluidThermo::Gamma() const
{
    tmp<volScalarField> tmpF
    (
        volScalarField::New
        (
            IOobject::groupName("Gamma", basicThermoModel::name_),
            volumeFractions_[0]/thermos_[0].Gamma()
        )
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]/thermos_[phasei].Gamma();
    }
    return 1.0/tmpF;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidThermo::Gamma(const label patchi) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       /thermos_[0].Gamma(patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           /thermos_[phasei].Gamma(patchi);
    }
    return 1.0/tmpF;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseFluidThermo::ESource() const
{
    tmp<volScalarField> tmpF
    (
        volScalarField::New
        (
            IOobject::groupName("ESource", basicThermoModel::name_),
            volumeFractions_[0]*thermos_[0].ESource()
        )
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].ESource();
    }
    return tmpF;
}



Foam::tmp<Foam::volScalarField> Foam::multiphaseFluidThermo::W() const
{
    tmp<volScalarField> tmpF
    (
        volScalarField::New
        (
            IOobject::groupName("W", basicThermoModel::name_),
            volumeFractions_[0]*thermos_[0].W()
        )
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].W();
    }
    return tmpF;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidThermo::W(const label patchi) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].W(patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].W(patchi);
    }
    return tmpF;
}


Foam::scalar Foam::multiphaseFluidThermo::Wi(const label celli) const
{
    scalar f(volumeFractions_[0][celli]*thermos_[0].Wi(celli));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        f += volumeFractions_[phasei][celli]*thermos_[phasei].Wi(celli);
    }
    return f;
}



Foam::tmp<Foam::volScalarField> Foam::multiphaseFluidThermo::Cp() const
{
    tmp<volScalarField> tmpF
    (
        volScalarField::New
        (
            IOobject::groupName("Cp", basicThermoModel::name_),
            volumeFractions_[0]*thermos_[0].Cp()
        )
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].Cp();
    }
    return tmpF;
}



Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidThermo::Cp(const label patchi) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].Cp(patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].Cp(patchi);
    }
    return tmpF;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidThermo::Cp
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].Cp(rho, e, T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].Cp(rho, e, T, patchi);
    }
    return tmpF;
}


Foam::scalar Foam::multiphaseFluidThermo::Cpi(const label celli) const
{
    scalar f(volumeFractions_[0][celli]*thermos_[0].Cpi(celli));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        f += volumeFractions_[phasei][celli]*thermos_[phasei].Cpi(celli);
    }
    return f;
}



Foam::tmp<Foam::volScalarField> Foam::multiphaseFluidThermo::Cv() const
{
    tmp<volScalarField> tmpF
    (
        volScalarField::New
        (
            IOobject::groupName("Cv", basicThermoModel::name_),
            volumeFractions_[0]*thermos_[0].Cv()
        )
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].Cv();
    }
    return tmpF;
}



Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidThermo::Cv(const label patchi) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].Cv(patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].Cv(patchi);
    }
    return tmpF;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidThermo::Cv
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].Cv(rho, e, T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].Cv(rho, e, T, patchi);
    }
    return tmpF;
}


Foam::scalar Foam::multiphaseFluidThermo::Cvi(const label celli) const
{
    scalar f(volumeFractions_[0][celli]*thermos_[0].Cvi(celli));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        f += volumeFractions_[phasei][celli]*thermos_[phasei].Cvi(celli);
    }
    return f;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseFluidThermo::CpByCv() const
{
    tmp<volScalarField> tmpF
    (
        volScalarField::New
        (
            IOobject::groupName("CpByCv", basicThermoModel::name_),
            volumeFractions_[0]*thermos_[0].CpByCv()
        )
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].CpByCv();
    }
    return tmpF;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidThermo::CpByCv(const label patchi) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].CpByCv(patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].CpByCv(patchi);
    }
    return tmpF;

}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidThermo::CpByCv
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].CpByCv(rho, e, T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].CpByCv(rho, e, T, patchi);
    }
    return tmpF;
}

// ************************************************************************* //
