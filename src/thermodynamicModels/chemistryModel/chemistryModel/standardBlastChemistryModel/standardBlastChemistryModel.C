/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "standardBlastChemistryModel.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::
standardBlastChemistryModel
(
    const blastThermo& thermo
)
:
    basicBlastChemistryModel(thermo),
    ODESystem(),
    mixture_
    (
        dynamic_cast<const mixtureBlastThermo<BasicThermo, ThermoType>&>(thermo)
    ),
    Y_(mixture_.composition().Y()),
    specieThermos_(mixture_.specieThermos()),
    reactions_(mixture_.species(), specieThermos_, this->mesh(), *this),
    nSpecie_(Y_.size()),
    nReaction_(reactions_.size()),
    Treact_(basicBlastChemistryModel::template lookupOrDefault<scalar>("Treact", 0)),
    RR_(nSpecie_),
    c_(nSpecie_),
    dcdt_(nSpecie_),
    cps_(nSpecie_),
    has_(nSpecie_)
{
    // Create the fields for the chemistry sources
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + Y_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mixture_.T().mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            )
        );
    }

    Info<< "standardBlastChemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::~standardBlastChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class ThermoType>
void Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcdt
) const
{
    dcdt = Zero;

    forAll(reactions_, i)
    {
        const Reaction<ThermoType>& R = reactions_[i];

        R.omega(p, T, c, li, dcdt);
    }
}


template<class BasicThermo, class ThermoType>
void Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::derivatives
(
    const scalar t,
    const scalarField& c,
    const label li,
    scalarField& dcdt
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    forAll(c_, i)
    {
        c_[i] = max(c[i], 0);
    }

    dcdt = Zero;

    // Evaluate contributions from reactions
    omega(p, T, c_, li, dcdt);

    // Evaluate the effect on the thermodynamic system ...

    // c*Cp & dT/dt
    scalar ccp = 0;
    scalar& dTdt = dcdt[nSpecie_];
    for (label i = 0; i < nSpecie_; i++)
    {
        const scalar rho = specieThermos_[i].rho(p, T);
        const scalar e = specieThermos_[i].Es(rho, T, T);

        ccp += c_[i]*specieThermos_[i].cp(rho, e, T);
        dTdt -= dcdt[i]*specieThermos_[i].ha(rho, e, T);
    }
    dTdt /= ccp;

    // dp/dt = 0 (pressure is assumed constant)
}


template<class BasicThermo, class ThermoType>
void Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    const label li,
    scalarField& dcdt,
    scalarSquareMatrix& J
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    forAll(c_, i)
    {
        c_[i] = max(c[i], 0);
    }

    dcdt = Zero;
    J = Zero;

    // Evaluate contributions from reactions
    forAll(reactions_, ri)
    {
        const Reaction<ThermoType>& R = reactions_[ri];
        scalar omegaI, kfwd, kbwd;
        const labelList null;
        R.dwdc(p, T, c_, li, J, dcdt, omegaI, kfwd, kbwd, false, null);
        R.dwdT(p, T, c_, li, omegaI, kfwd, kbwd, J, false, null, nSpecie_);
    }

    // Evaluate the effect on the thermodynamic system ...

    // c*Cp & dT/dt
    scalar ccp = 0, dccpdT = 0;
    scalar& dTdt = dcdt[nSpecie_];
    for (label i = 0; i < nSpecie_; i++)
    {
        const scalar rho = specieThermos_[i].rho(p, T);
        const scalar e = specieThermos_[i].Es(rho, T, T);

        cps_[i] = specieThermos_[i].cp(rho, e, T);
        has_[i] = specieThermos_[i].ha(rho, e, T);

        ccp += c_[i]*cps_[i];
        dccpdT += c_[i]*specieThermos_[i].dcpdT(rho, e, T);
        dTdt -= dcdt[i]*has_[i];
    }
    dTdt /= ccp;


    // dp/dt = 0 (pressure is assumed constant)

    // d(dTdt)/dc
    for (label i = 0; i < nSpecie_; i++)
    {
        scalar& d2Tdtdci = J(nSpecie_, i);
        for (label j = 0; j < nSpecie_; j++)
        {
            const scalar d2cjdtdci = J(j, i);
            d2Tdtdci -= d2cjdtdci*has_[j];
        }
        d2Tdtdci -= cps_[i]*dTdt;
        d2Tdtdci /= ccp;
    }

    // d(dTdt)/dT
    scalar& d2TdtdT = J(nSpecie_, nSpecie_);
    for (label i = 0; i < nSpecie_; i++)
    {
        const scalar d2cidtdT = J(i, nSpecie_);
        d2TdtdT -= dcdt[i]*cps_[i] + d2cidtdT*has_[i];
    }
    d2TdtdT -= dTdt*dccpdT;
    d2TdtdT /= ccp;

    // d(dpdt)/dc = 0 (pressure is assumed constant)

    // d(dpdt)/dT = 0 (pressure is assumed constant)
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::tc() const
{
    tmp<volScalarField> ttc
    (
        volScalarField::New
        (
            "tc",
            this->mesh(),
            dimensionedScalar(dimTime, small),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    scalarField& tc = ttc.ref();

    const scalarField& T = this->thermo().T();
    tmp<volScalarField> trho = this->thermo().rho();
    const scalarField& rho = trho();

    if (this->chemistry_)
    {
        reactionEvaluationScope scope(*this);

        forAll(rho, celli)
        {
            const scalar rhoi = rho[celli];
            const scalar Ti = T[celli];
            const scalar pi = this->mixture().cellp(celli);

            for (label i=0; i<nSpecie_; i++)
            {
                c_[i] = rhoi*Y_[i][celli]/specieThermos_[i].W();
            }

            // A reaction's rate scale is calculated as it's molar
            // production rate divided by the total number of moles in the
            // system.
            //
            // The system rate scale is the average of the reactions' rate
            // scales weighted by the reactions' molar production rates. This
            // weighting ensures that dominant reactions provide the largest
            // contribution to the system rate scale.
            //
            // The system time scale is then the reciprocal of the system rate
            // scale.
            //
            // Contributions from forward and reverse reaction rates are
            // handled independently and identically so that reversible
            // reactions produce the same result as the equivalent pair of
            // irreversible reactions.

            scalar sumW = 0, sumWRateByCTot = 0;
            forAll(reactions_, i)
            {
                const Reaction<ThermoType>& R = reactions_[i];
                scalar omegaf, omegar;
                R.omega(pi, Ti, c_, celli, omegaf, omegar);

                scalar wf = 0;
                forAll(R.rhs(), s)
                {
                    wf += R.rhs()[s].stoichCoeff*omegaf;
                }
                sumW += wf;
                sumWRateByCTot += sqr(wf);

                scalar wr = 0;
                forAll(R.lhs(), s)
                {
                    wr += R.lhs()[s].stoichCoeff*omegar;
                }
                sumW += wr;
                sumWRateByCTot += sqr(wr);
            }

            tc[celli] =
                sumWRateByCTot == 0 ? vGreat : sumW/sumWRateByCTot*sum(c_);
        }
    }

    ttc.ref().correctBoundaryConditions();

    return ttc;
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        volScalarField::New
        (
            "Qdot",
            this->mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        reactionEvaluationScope scope(*this);

        scalarField& Qdot = tQdot.ref();

        forAll(Y_, i)
        {
            forAll(Qdot, celli)
            {
                const scalar hi = specieThermos_[i].Hf();
                Qdot[celli] -= hi*RR_[i][celli];
            }
        }
    }

    return tQdot;
}


template<class BasicThermo, class ThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::calculateRR
(
    const label ri,
    const label si
) const
{
    tmp<volScalarField::Internal> tRR
    (
        volScalarField::Internal::New
        (
            "RR",
            this->mesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );

    volScalarField::Internal& RR = tRR.ref();

    const scalarField& T = this->thermo().T();
    tmp<volScalarField> trho = this->thermo().rho();
    const scalarField& rho = trho();

    reactionEvaluationScope scope(*this);

    scalar omegaf, omegar;

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = this->mixture().cellp(celli);

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c_[i] = rhoi*Yi/specieThermos_[i].W();
        }

        const Reaction<ThermoType>& R = reactions_[ri];
        const scalar omegaI = R.omega(pi, Ti, c_, celli, omegaf, omegar);

        forAll(R.lhs(), s)
        {
            if (si == R.lhs()[s].index)
            {
                RR[celli] -= R.lhs()[s].stoichCoeff*omegaI;
            }
        }

        forAll(R.rhs(), s)
        {
            if (si == R.rhs()[s].index)
            {
                RR[celli] += R.rhs()[s].stoichCoeff*omegaI;
            }
        }

        RR[celli] *= specieThermos_[si].W();
    }

    return tRR;
}


template<class BasicThermo, class ThermoType>
void Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    const scalarField& T = this->thermo().T();
    tmp<volScalarField> trho = this->thermo().rho();
    const scalarField& rho = trho();

    reactionEvaluationScope scope(*this);

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = this->mixture().cellp(celli);

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c_[i] = rhoi*Yi/specieThermos_[i].W();
        }

        omega(pi, Ti, c_, celli, dcdt_);

        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] = dcdt_[i]*specieThermos_[i].W();
        }
    }
}


template<class BasicThermo, class ThermoType>
template<class DeltaTType>
Foam::scalar
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    basicBlastChemistryModel::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    const scalarField& T = this->thermo().T();
    tmp<volScalarField> trho0 = this->thermo().rho();
    const scalarField& rho0 = trho0();

    reactionEvaluationScope scope(*this);

    scalarField c0(nSpecie_);

    forAll(rho0, celli)
    {
        scalar Ti = T[celli];

        if (Ti > Treact_)
        {
            const scalar rhoi = rho0[celli];
            scalar pi = this->mixture().cellp(celli);

            for (label i=0; i<nSpecie_; i++)
            {
                c_[i] = rhoi*Y_[i].oldTime()[celli]/specieThermos_[i].W();
                c0[i] = c_[i];
            }

            // Initialise time progress
            scalar timeLeft = deltaT[celli];

            // Calculate the chemical source terms
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                this->solve(pi, Ti, c_, celli, dt, this->deltaTChem_[celli]);
                timeLeft -= dt;
            }

            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            this->deltaTChem_[celli] =
                min(this->deltaTChem_[celli], this->deltaTChemMax_);

            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] =
                    (c_[i] - c0[i])*specieThermos_[i].W()/deltaT[celli];
            }
        }
        else
        {
            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = 0;
            }
        }
    }

    return deltaTMin;
}


template<class BasicThermo, class ThermoType>
Foam::scalar Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class BasicThermo, class ThermoType>
Foam::scalar
Foam::standardBlastChemistryModel<BasicThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


// ************************************************************************* //
