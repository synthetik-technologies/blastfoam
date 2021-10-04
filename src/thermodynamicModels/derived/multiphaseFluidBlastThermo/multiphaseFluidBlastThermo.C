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

#include "multiphaseFluidBlastThermo.H"
#include "scalarEquation.H"
#include "NewtonRaphsonRootSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphaseFluidBlastThermo, 0);
    addToRunTimeSelectionTable
    (
        fluidBlastThermo,
        multiphaseFluidBlastThermo,
        phase
    );
}


// * * * * * * * * * multiphaseEEquation Member Functions  * * * * * * * * * //

namespace Foam
{
    class multiphaseEEquation
    :
        public scalarEquation
    {
        const multiphaseFluidBlastThermo& thermo_;
        const volScalarField& p_;
        mutable volScalarField* e_;

    public:
        multiphaseEEquation
        (
            multiphaseFluidBlastThermo& thermo
        )
        :
            thermo_(thermo),
            p_(thermo_.p()),
            e_(&thermo.he())
        {}
        virtual scalar f(const scalar e, const label li) const
        {
            (*e_)[li] = e;
            return thermo_.cellpRhoT(li, false) - p_[li];
        }
        virtual scalar dfdx(const scalar e, const label li) const
        {
            (*e_)[li] = e;
            return thermo_.celldpde(li);
        }
    };

    class multiphaseTHEEquation
    :
        public scalarEquation
    {
        const multiphaseFluidBlastThermo& thermo_;
        const volScalarField& he_;
        label patchi_;

    public:
        multiphaseTHEEquation
        (
            multiphaseFluidBlastThermo& thermo,
            const scalar TLow
        )
        :
            scalarEquation(TLow, great),
            thermo_(thermo),
            he_(thermo.he()),
            patchi_(-1)
        {}

        label& patch()
        {
            return patchi_;
        }

        virtual scalar f(const scalar T, const label li) const
        {
            return
                patchi_ == -1
              ? thermo_.cellHE(T, li) - he_[li]
              : thermo_.patchFaceHE(T, patchi_, li)
              - he_.boundaryField()[patchi_][li];
        }
        virtual scalar dfdx(const scalar T, const label li) const
        {
            return
                patchi_ == -1
              ? thermo_.cellCpv(T, li)
              : thermo_.patchFaceCpv(T, patchi_, li);
        }
    };
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiphaseFluidBlastThermo::calculate()
{
    scalarField& TCells = T_.primitiveFieldRef();
    scalarField& heCells = this->he().primitiveFieldRef();
    volScalarField::Boundary& bT = T_.boundaryFieldRef();
    volScalarField::Boundary& bhe = this->he().boundaryFieldRef();

    multiphaseTHEEquation eqn(*this, this->TLow_);
    NewtonRaphsonRootSolver solver(eqn, dictionary());
    forAll(TCells, celli)
    {
        TCells[celli] = solver.solve(TCells[celli], celli);
    }
    forAll(bT, patchi)
    {
        eqn.patch() = patchi;
        scalarField& pT = bT[patchi];
        forAll(pT, facei)
        {
            pT[facei] = solver.solve(pT[facei], facei);
        }
    }
    T_.correctBoundaryConditions();

    if (min(T_).value() < this->TLow_)
    {
        forAll(TCells, celli)
        {
            if (TCells[celli] < this->TLow_)
            {
                TCells[celli] = this->TLow_;
                heCells[celli] = this->cellHE(this->TLow_, celli);
            }
        }
        forAll(bhe, patchi)
        {
            scalarField& pT = bT[patchi];
            scalarField& phe = bhe[patchi];

            pT = max(pT, this->TLow_);
            phe = this->he(pT, patchi);
        }
    }
    this->he().correctBoundaryConditions();

    volScalarField XiSum
    (
        volScalarField::New("XiSum", mesh(), dimensionedScalar(dimless, 0.0))
    );
    volScalarField pXiSum
    (
        volScalarField::New
        (
            "pXiSum",
            mesh(),
            dimensionedScalar(dimPressure, 0.0)
        )
    );
    volScalarField cSqrRhoXiSum
    (
        volScalarField::New
        (
            "cSqrRhoXiSum",
            mesh(),
            dimensionedScalar(sqr(dimVelocity)*dimDensity, 0.0)
        )
    );

    this->Cp_ = Zero;
    this->Cv_ = Zero;
    this->alpha_ = Zero;
    this->mu_ = Zero;

    forAll(thermos_, i)
    {
        thermos_[i].calculate
        (
            volumeFractions_[i],
            this->he(),
            this->T_,
            this->Cp_,
            this->Cv_,
            this->alpha_,
            this->mu_,
            pXiSum,
            XiSum
        );
    }
    XiSum.max(small);
    this->p_ = pXiSum/XiSum;
    this->p_.correctBoundaryConditions();

    // Normalise variables
    normalise(Cp_);
    normalise(Cv_);
    normalise(mu_);
    normalise(alpha_);

    forAll(thermos_, i)
    {
        thermos_[i].calculateSpeedOfSound
        (
            volumeFractions_[i],
            cSqrRhoXiSum
        );
    }
    cSqrRhoXiSum.max(small);
    this->speedOfSound_ =
        sqrt
        (
            cSqrRhoXiSum
           /max(this->rho_*XiSum, this->residualRho_)
        );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseFluidBlastThermo::multiphaseFluidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    fluidBlastThermo(mesh, dict, phaseName),
    phases_(dict.lookup("phases")),
    volumeFractions_(phases_.size()),
    rhos_(phases_.size()),
    thermos_(phases_.size()),
    alphaRhos_(phases_.size()),
    residualAlpha_(dimless, 0.0),
    sumVfPtr_(nullptr)

{
    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        0.0
    );
    forAll(phases_, phasei)
    {
        word phaseIName = phases_[phasei];
        volumeFractions_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alpha", phases_[phasei]),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
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
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
        thermos_.set
        (
            phasei,
            phaseFluidBlastThermo::New
            (
                mesh,
                dict.subDict(phaseIName),
                phaseIName,
                phaseName
            ).ptr()
        );

        thermos_[phasei].read(dict.subDict(phaseIName));
        residualAlpha_ = max(thermos_[phasei].residualAlpha(), residualAlpha_);
        residualRho_ = max(thermos_[phasei].residualRho(), residualRho_);

        rho_ += volumeFractions_[phasei]*rhos_[phasei];
        sumAlpha += volumeFractions_[phasei];
    }
    rho_ /= max(sumAlpha, residualAlpha_);
    this->initializeFields();
}


void Foam::multiphaseFluidBlastThermo::initializeModels()
{
    forAll(thermos_, phasei)
    {
        thermos_[phasei].initializeModels();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseFluidBlastThermo::~multiphaseFluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiphaseFluidBlastThermo::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


void Foam::multiphaseFluidBlastThermo::read(const dictionary& dict)
{
    forAll(thermos_, phasei)
    {
        thermos_[phasei].read(dict.subDict(thermos_[phasei].phaseName()));
        residualAlpha_ = max(residualAlpha_, thermos_[phasei].residualAlpha());
    }
}


void Foam::multiphaseFluidBlastThermo::setTotalVolumeFractionPtr
(
    const volScalarField& vf
)
{
    sumVfPtr_.set(&vf);
}


void Foam::multiphaseFluidBlastThermo::update()
{
    forAll(phases_, phasei)
    {
        thermos_[phasei].update();
    }
}


void Foam::multiphaseFluidBlastThermo::solve()
{
    forAll(phases_, phasei)
    {
        thermos_[phasei].solve();
    }
}


void Foam::multiphaseFluidBlastThermo::postUpdate()
{
    forAll(phases_, phasei)
    {
        thermos_[phasei].postUpdate();
    }
}


void Foam::multiphaseFluidBlastThermo::updateRho(const volScalarField& p)
{
    thermos_[0].updateRho(p);
    this->rho_ = volumeFractions_[0]*thermos_[0].rho();
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        thermos_[phasei].updateRho(p);
        rho_ += volumeFractions_[phasei]*thermos_[phasei].rho();
    }
    normalise(rho_);
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseFluidBlastThermo::ESource() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].ESource());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].ESource();
    }
    return tmpF;
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::calce(const volScalarField& p) const
{
    tmp<volScalarField> eInitTmp(volScalarField::New("eInit", e_));
    volScalarField& eInit(eInitTmp.ref());
    multiphaseFluidBlastThermo& thermo
    (
        const_cast<multiphaseFluidBlastThermo&>(*this)
    );
    volScalarField& e(const_cast<volScalarField&>(e_));

    multiphaseEEquation eqn(thermo);
    dictionary dict;
    dict.add("tolerance", 1e-6);
    NewtonRaphsonRootSolver solver(eqn, dict);
    forAll(eInit, celli)
    {
        if (mag(celldpde(celli)) < small)
        {
            eInit[celli] = cellHE(T_[celli], celli);
        }
        else
        {
            const scalar e0 = e[celli];
            eInit[celli] = solver.solve(e0, celli);
            e[celli] = e0;
        }
    }

    return eInitTmp;
}


Foam::scalar Foam::multiphaseFluidBlastThermo::cellpRhoT
(
    const label celli,
    const bool limit
) const
{
    scalar rGamma = 0.0;
    scalar pByGamma = 0.0;
    forAll(thermos_, phasei)
    {
        scalar alphai(volumeFractions_[phasei][celli]);
        if (alphai > residualAlpha_.value())
        {
            scalar Xi(alphai/(thermos_[phasei].cellGamma(celli) - 1.0));
            rGamma += Xi;
            pByGamma += Xi*thermos_[phasei].cellpRhoT(celli, limit);
        }
    }
    return pByGamma/max(rGamma, residualAlpha_.value());
}


Foam::scalar Foam::multiphaseFluidBlastThermo::cellGamma(const label celli) const
{
    scalar f(volumeFractions_[0][celli]/thermos_[0].cellGamma(celli));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        f += volumeFractions_[phasei][celli]/thermos_[phasei].cellGamma(celli);
    }
    return 1.0/cellNormalise(f, celli);
}


Foam::scalar Foam::multiphaseFluidBlastThermo::celldpdRho(const label celli) const
{
    scalar rGamma = 0.0;
    scalar dpdRho = 0.0;
    forAll(thermos_, phasei)
    {
        scalar alphai(volumeFractions_[phasei][celli]);
        if (alphai > residualAlpha_.value())
        {
            scalar Xi(alphai/(thermos_[phasei].cellGamma(celli) - 1.0));
            rGamma += Xi;
            dpdRho += Xi*thermos_[phasei].celldpdRho(celli);
        }
    }
    return dpdRho/max(rGamma, residualAlpha_.value());
}


Foam::scalar Foam::multiphaseFluidBlastThermo::celldpde(const label celli) const
{
    scalar rGamma = 0.0;
    scalar dpde = 0.0;
    forAll(thermos_, phasei)
    {
        scalar alphai(volumeFractions_[phasei][celli]);
        if (alphai > residualAlpha_.value())
        {
            scalar Xi(alphai/(thermos_[phasei].cellGamma(celli) - 1.0));
            rGamma += Xi;
            dpde += Xi*thermos_[phasei].celldpde(celli);
        }
    }
    return dpde/max(rGamma, residualAlpha_.value());
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].he(p, T));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].he(p, T);
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> tmpF
    (
        UIndirectList<scalar>(volumeFractions_[0], cells)()
       *thermos_[0].he(T, cells)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            UIndirectList<scalar>(volumeFractions_[phasei], cells)()
           *thermos_[phasei].he(T, cells);
    }
    return normalise(tmpF, cells);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::he
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].he(T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].he(T, patchi);
    }
    return normalise(tmpF, patchi);
}


Foam::scalar Foam::multiphaseFluidBlastThermo::cellHE
(
    const scalar T,
    const label celli
) const
{
    scalar heSum = 0.0;
    forAll(thermos_, phasei)
    {
        scalar alphai(volumeFractions_[phasei][celli]);
        if (alphai > thermos_[phasei].residualAlpha().value())
        {
            heSum += alphai*thermos_[phasei].cellHE(T, celli);
        }
    }
    return cellNormalise(heSum, celli);
}


Foam::scalar Foam::multiphaseFluidBlastThermo::patchFaceHE
(
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalar heSum = 0.0;
    forAll(thermos_, phasei)
    {
        scalar alphai
        (
            volumeFractions_[phasei].boundaryField()[patchi][facei]
        );
        if (alphai > thermos_[phasei].residualAlpha().value())
        {
            heSum += alphai*thermos_[phasei].patchFaceHE(T, patchi, facei);
        }
    }
    return patchFaceNormalise(heSum, patchi, facei);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::hs() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].hs());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].hs();
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].hs(p, T));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].hs(p, T);
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> tmpF
    (
        UIndirectList<scalar>(volumeFractions_[0], cells)()
       *thermos_[0].hs(T, cells)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            UIndirectList<scalar>(volumeFractions_[phasei], cells)()
           *thermos_[phasei].hs(T, cells);
    }
    return normalise(tmpF, cells);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].hs(T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].hs(T, patchi);
    }
    return normalise(tmpF, patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::ha() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].ha());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].ha();
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].ha(p, T));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].ha(p, T);
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> tmpF
    (
        UIndirectList<scalar>(volumeFractions_[0], cells)()
       *thermos_[0].ha(T, cells)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            UIndirectList<scalar>(volumeFractions_[phasei], cells)()
           *thermos_[phasei].ha(T, cells);
    }
    return normalise(tmpF, cells);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::ha
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].ha(T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].ha(T, patchi);
    }
    return normalise(tmpF, patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::hc() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].hc());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].hc();
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::flameT() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].flameT());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].flameT();
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::THE() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].THE());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].THE();
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::THE
(
    const volScalarField& he,
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].THE(he, p, T));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].THE(he, p, T);
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::THE
(
    const scalarField& he,
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> tmpF
    (
        UIndirectList<scalar>(volumeFractions_[0], cells)()
       *thermos_[0].THE(he, T, cells)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            UIndirectList<scalar>(volumeFractions_[phasei], cells)()
           *thermos_[phasei].THE(he, T, cells);
    }
    return normalise(tmpF, cells);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::THE
(
    const scalarField& he,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].THE(he, T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].THE(he, T, patchi);
    }
    return normalise(tmpF, patchi);
}


Foam::scalar Foam::multiphaseFluidBlastThermo::cellTHE
(
    const scalar he,
    const scalar T0,
    const label celli
) const
{
    scalar f(volumeFractions_[0][celli]*thermos_[0].cellTHE(he, T0, celli));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        f +=
            volumeFractions_[phasei][celli]
           *thermos_[phasei].cellTHE(he, T0, celli);
    }
    return cellNormalise(f, celli);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].Cp(T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].Cp(T, patchi);
    }
    return normalise(tmpF, patchi);
}


Foam::scalar Foam::multiphaseFluidBlastThermo::cellCp
(
    const scalar T,
    const label celli
) const
{
    scalar f(volumeFractions_[0][celli]*thermos_[0].cellCp(T, celli));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        f += volumeFractions_[phasei][celli]*thermos_[phasei].cellCp(T, celli);
    }
    return cellNormalise(f, celli);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].Cp(T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].Cv(T, patchi);
    }
    return normalise(tmpF, patchi);
}


Foam::scalar Foam::multiphaseFluidBlastThermo::cellCv
(
    const scalar T,
    const label celli
) const
{
    scalar f(volumeFractions_[0][celli]*thermos_[0].cellCv(T, celli));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        f += volumeFractions_[phasei][celli]*thermos_[phasei].cellCv(T, celli);
    }
    return cellNormalise(f, celli);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].Cpv(T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].Cpv(T, patchi);
    }
    return normalise(tmpF, patchi);
}


Foam::scalar Foam::multiphaseFluidBlastThermo::cellCpv
(
    const scalar T,
    const label celli
) const
{
    scalar cpvSum = 0.0;
    forAll(thermos_, phasei)
    {
        scalar alphai(volumeFractions_[phasei][celli]);
        if (alphai > residualAlpha_.value())
        {
            cpvSum += alphai*thermos_[phasei].cellCpv(T, celli);
        }
    }
    return cellNormalise(cpvSum, celli);
}


Foam::scalar Foam::multiphaseFluidBlastThermo::patchFaceCpv
(
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalar cpvSum = 0.0;
    forAll(thermos_, phasei)
    {
        scalar alphai
        (
            volumeFractions_[phasei].boundaryField()[patchi][facei]
        );
        if (alphai > residualAlpha_.value())
        {
            cpvSum += alphai*thermos_[phasei].patchFaceCpv(T, patchi, facei);
        }
    }
    return patchFaceNormalise(cpvSum, patchi, facei);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::W() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].W());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].W();
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::W(const label patchi) const
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
    return 1.0/normalise(tmpF, patchi);
}


Foam::scalar Foam::multiphaseFluidBlastThermo::cellW
(
    const label celli
) const
{
    scalar f(volumeFractions_[0][celli]*thermos_[0].cellW(celli));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        f += volumeFractions_[phasei][celli]*thermos_[phasei].cellW(celli);
    }
    return cellNormalise(f, celli);
}


// ************************************************************************* //
