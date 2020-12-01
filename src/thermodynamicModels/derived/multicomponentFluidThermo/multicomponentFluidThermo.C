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

#include "multicomponentFluidThermo.H"
#include "fvc.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
const Thermo& Foam::multicomponentFluidThermo<Thermo>::constructSpeciesData
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new Thermo(thermoDict.subDict(species_[i]))
        );
    }

    return speciesData_[0];
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multicomponentFluidThermo<Thermo>::multicomponentFluidThermo
(
    const word& name,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master,
    const word& masterName
)
:
    fluidThermoModel
    (
        name,
        p,
        rho,
        e,
        T,
        dict,
        master,
        masterName
    ),
    species_(dict.lookup("species")),
    Ys_(species_.size()),
    speciesData_(species_.size()),
    mixture_(constructSpeciesData(dict)),
    YsOld_(species_.size()),
    deltaAlphaRhoYs_(species_.size())
{
    tmp<volScalarField> tYdefault;
    volScalarField YTot(volScalarField::New("YTot", p.mesh(), 0.0));

    const fvMesh& mesh = rho.mesh();
    forAll(species_, i)
    {
        IOobject header
        (
            IOobject::groupName(species_[i], name),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        // check if field exists and can be read
        if (header.typeHeaderOk<volScalarField>(true))
        {
            Ys_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(species_[i], name),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );
        }
        else
        {
            // Read Ydefault if not already read
            if (!tYdefault.valid())
            {
                word YdefaultName(IOobject::groupName("Ydefault", name));

                IOobject timeIO
                (
                    YdefaultName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                IOobject time0IO
                (
                    YdefaultName,
                    Time::timeName(0),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                if (timeIO.typeHeaderOk<volScalarField>(true))
                {
                    tYdefault = new volScalarField(timeIO, mesh);
                }
                else
                {
                    tYdefault = new volScalarField(time0IO, mesh);
                }
            }

            Ys_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName(species_[i], name),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    tYdefault()
                )
            );
        }
        YTot += Ys_[i];
    }
    forAll(Ys_, i)
    {
        Ys_[i] /= YTot;
        Ys_[i].correctBoundaryConditions();
    }

    updateMixture();

    //- Initialize the density using the pressure and temperature
    //  This is only done at the first time step (Not on restart)
    if
    (
        max(this->rho_).value() == 0
     || (
            dict.lookupOrDefault<Switch>("calculateDensity", false)
         && this->rho_.time().timeIndex() == 0
        )
    )
    {
        forAll(this->T_, celli)
        {
            const Thermo& mixture(volMixture_[celli]);
            this->rho_[celli] = mixture.initializeRho
            (
                this->p_[celli],
                this->rho_[celli],
                this->e_[celli],
                this->T_[celli]
            );
        }
        forAll(this->T_.boundaryField(), patchi)
        {
            forAll(this->T_.boundaryField()[patchi], facei)
            {
                const Thermo& mixture(faceMixture_[patchi][facei]);
                this->rho_.boundaryFieldRef()[patchi][facei] =
                    mixture.initializeRho
                    (
                        this->p_.boundaryField()[patchi][facei],
                        this->rho_.boundaryField()[patchi][facei],
                        this->e_.boundaryField()[patchi][facei],
                        this->T_.boundaryField()[patchi][facei]
                    );
            }
        }
    }

    forAll(this->T_, celli)
        {
            const Thermo& mixture(volMixture_[celli]);
            this->fluidThermoModel::mu_[celli] = mixture.mu
            (
                this->rho_[celli],
                this->e_[celli],
                this->T_[celli]
            );
        }
        forAll(this->T_.boundaryField(), patchi)
        {
            forAll(this->T_.boundaryField()[patchi], facei)
            {
                const Thermo& mixture(faceMixture_[patchi][facei]);
                this->fluidThermoModel::mu_.boundaryFieldRef()[patchi][facei] =
                    mixture.mu
                    (
                        this->rho_.boundaryField()[patchi][facei],
                        this->e_.boundaryField()[patchi][facei],
                        this->T_.boundaryField()[patchi][facei]
                    );
            }
        }

    this->initialize();

    forAll(YsOld_, i)
    {
        YsOld_.set(i, new PtrList<volScalarField>());
        deltaAlphaRhoYs_.set(i, new PtrList<volScalarField>());
    }
    this->lookupAndInitialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multicomponentFluidThermo<Thermo>::~multicomponentFluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::multicomponentFluidThermo<Thermo>::initializeModels()
{
    if
    (
        this->rho_.mesh().template foundObject<surfaceScalarField>
        (
            IOobject::groupName("alphaRhoPhi", this->name_)
        )
    )
    {
        alphaRhoPhiPtr_ =
            &this->rho_.mesh().template lookupObject<surfaceScalarField>
            (
                IOobject::groupName("alphaRhoPhi", this->name_)
            );
    }
    else if
    (
        this->rho_.mesh().template foundObject<surfaceScalarField>
        (
            IOobject::groupName("alphaRhoPhi", this->masterName_)
        )
    )
    {
        alphaRhoPhiPtr_ =
            &this->rho_.mesh().template lookupObject<surfaceScalarField>
            (
                IOobject::groupName("alphaRhoPhi", this->masterName_)
            );
    }
    else
    {
        alphaRhoPhiPtr_ =
            &this->rho_.mesh().template lookupObject<surfaceScalarField>
            (
                IOobject::groupName("rhoPhi", this->masterName_)
            );
    }

    if
    (
        this->rho_.mesh().template foundObject<volScalarField>
        (
            IOobject::groupName("alphaRho", this->name_)
        )
    )
    {
        alphaRhoPtr_ =
            &this->rho_.mesh().template lookupObject<volScalarField>
            (
                IOobject::groupName("alphaRho", this->name_)
            );
    }
    else if
    (
        this->rho_.mesh().template foundObject<volScalarField>
        (
            IOobject::groupName("alphaRho", this->masterName_)
        )
    )
    {
        alphaRhoPtr_ =
            &this->rho_.mesh().template lookupObject<volScalarField>
            (
                IOobject::groupName("alphaRho", this->masterName_)
            );
    }
    else
    {
        alphaRhoPtr_ =
            &this->rho_.mesh().template lookupObject<volScalarField>
            (
                IOobject::groupName("rho", this->masterName_)
            );
    }
}


template<class Thermo>
void Foam::multicomponentFluidThermo<Thermo>::solve()
{
    const dimensionedScalar& dT(this->rho_.mesh().time().deltaT());
    volScalarField YTot(volScalarField::New("YTot", this->rho_.mesh(), 0.0));
    forAll(Ys_, i)
    {
        volScalarField YOld(Ys_[i]);
        this->storeAndBlendOld(YOld, YsOld_[i]);

        volScalarField deltaAlphaRhoY
        (
            fvc::div(*alphaRhoPhiPtr_, Ys_[i], "div(Y)")
        );
        this->storeAndBlendDelta(deltaAlphaRhoY, deltaAlphaRhoYs_[i]);

        Ys_[i] =
            (
                (*alphaRhoPtr_).prevIter()*YOld - dT*deltaAlphaRhoY
            )/max(this->residualRho(), *alphaRhoPtr_);
        Ys_[i].maxMin(0.0, 1.0);

        YTot += Ys_[i];
    }
    forAll(Ys_, i)
    {
        Ys_[i] /= YTot;
        Ys_[i].correctBoundaryConditions();
    }
}


template<class Thermo>
void Foam::multicomponentFluidThermo<Thermo>::clearODEFields()
{
    forAll(YsOld_, i)
    {
        this->clearOld(YsOld_[i]);
        this->clearDelta(deltaAlphaRhoYs_[i]);
    }
}



template<class Thermo>
void Foam::multicomponentFluidThermo<Thermo>::updateMixture()
{
    volMixture_.setSize(this->rho_.size());
    forAll(volMixture_, celli)
    {
        volMixture_.set(celli, new Thermo(cellMixture(celli)));
    }

    faceMixture_.setSize(this->rho_.boundaryField().size());
    forAll(faceMixture_, patchi)
    {
        faceMixture_.set
        (
            patchi,
            new PtrList<Thermo>(this->rho_.boundaryField()[patchi].size())
        );
        forAll(faceMixture_[patchi], facei)
        {
            faceMixture_[patchi].set
            (
                facei,
                new Thermo(patchFaceMixture(patchi, facei))
            );
        }
    }
}


template<class Thermo>
const Thermo& Foam::multicomponentFluidThermo<Thermo>::cellMixture
(
    const label celli
) const
{
    mixture_ = Ys_[0][celli]*speciesData_[0];

    for (label n = 1; n < Ys_.size(); n++)
    {
        mixture_ += Ys_[n][celli]*speciesData_[n];
    }

    return mixture_;
}


template<class Thermo>
const Thermo& Foam::multicomponentFluidThermo<Thermo>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    mixture_ = Ys_[0].boundaryField()[patchi][facei]*speciesData_[0];

    for (label n = 1; n < Ys_.size(); n++)
    {
        mixture_ += Ys_[n].boundaryField()[patchi][facei]*speciesData_[n];
    }

    return mixture_;
}


template<class Thermo>
void Foam::multicomponentFluidThermo<Thermo>::correct()
{
    updateMixture();
    forAll(this->T_, celli)
    {
        const Thermo& mixture(volMixture_[celli]);
        if (this->master_)
        {
            this->T_[celli] = mixture.TRhoE
            (
                this->T_[celli],
                this->rho_[celli],
                this->e_[celli]
            );
            this->p_[celli] = mixture.p
            (
                this->rho_[celli],
                this->e_[celli],
                this->T_[celli]
            );
        }

        if (this->viscous_)
        {
            this->fluidThermoModel::mu_[celli] = mixture.mu
            (
                this->rho_[celli],
                this->e_[celli],
                this->T_[celli]
            );
            this->fluidThermoModel::alpha_[celli] = mixture.alphah
            (
                this->rho_[celli],
                this->e_[celli],
                this->T_[celli]
            );
        }
    }
    forAll(this->T_.boundaryField(), patchi)
    {
        forAll(this->T_.boundaryField()[patchi], facei)
        {
            const Thermo& mixture(faceMixture_[patchi][facei]);
            if (this->master_)
            {
                this->T_.boundaryFieldRef()[patchi][facei] = mixture.TRhoE
                (
                    this->T_.boundaryField()[patchi][facei],
                    this->rho_.boundaryField()[patchi][facei],
                    this->e_.boundaryField()[patchi][facei]
                );
                this->p_.boundaryFieldRef()[patchi][facei] = mixture.p
                (
                    this->rho_.boundaryField()[patchi][facei],
                    this->e_.boundaryField()[patchi][facei],
                    this->T_.boundaryField()[patchi][facei]
                );
            }

            if (this->viscous_)
            {
                this->fluidThermoModel::mu_.boundaryFieldRef()[patchi][facei] =
                mixture.mu
                    (
                        this->rho_.boundaryField()[patchi][facei],
                        this->e_.boundaryField()[patchi][facei],
                        this->T_.boundaryField()[patchi][facei]
                    );
                this->fluidThermoModel::alpha_.boundaryFieldRef()[patchi][facei] =
                    mixture.alphah
                    (
                        this->rho_.boundaryField()[patchi][facei],
                        this->e_.boundaryField()[patchi][facei],
                        this->T_.boundaryField()[patchi][facei]
                    );
            }
        }
    }
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidThermo<Thermo>::ESource() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ESource",
                this->p_.mesh().time().timeName(),
                this->p_.mesh()
            ),
            this->p_.mesh(),
            dimensionedScalar("0", dimEnergy/dimTime/dimVolume, 0.0)
        )
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidThermo<Thermo>::speedOfSound() const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName("speedOfSound", this->group()),
            this->p_.mesh(),
            dimVelocity
        )
    );
    volScalarField& psi(tPsi.ref());
    forAll(this->T_, celli)
    {
        const Thermo& mixture(volMixture_[celli]);
        psi[celli] = mixture.speedOfSound
        (
            this->p_[celli],
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        );
    }
    forAll(this->T_.boundaryField(), patchi)
    {
        forAll(this->T_.boundaryField()[patchi], facei)
        {
            const Thermo& mixture(faceMixture_[patchi][facei]);
            psi.boundaryFieldRef()[patchi][facei] =
                mixture.speedOfSound
                (
                    this->p_.boundaryField()[patchi][facei],
                    this->rho_.boundaryField()[patchi][facei],
                    this->e_.boundaryField()[patchi][facei],
                    this->T_.boundaryField()[patchi][facei]
                );
        }
    }
    return tPsi;
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidThermo<Thermo>::speedOfSound(const label patchi) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->T_.boundaryField()[patchi], facei)
    {
        const Thermo& mixture(faceMixture_[patchi][facei]);
        psi[facei] =
            mixture.speedOfSound
            (
                this->p_.boundaryField()[patchi][facei],
                this->rho_.boundaryField()[patchi][facei],
                this->e_.boundaryField()[patchi][facei],
                this->T_.boundaryField()[patchi][facei]
            );
    }
    return tPsi;
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidThermo<Thermo>::calcP(const label patchi) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->T_.boundaryField()[patchi], facei)
    {
        const Thermo& mixture(faceMixture_[patchi][facei]);
        psi[facei] =
            mixture.p
            (
                this->rho_.boundaryField()[patchi][facei],
                this->e_.boundaryField()[patchi][facei],
                this->T_.boundaryField()[patchi][facei]
            );
    }
    return tPsi;
}


template<class Thermo>
Foam::scalar Foam::multicomponentFluidThermo<Thermo>::calcPi(const label celli) const
{
    return
        volMixture_[celli].p
        (
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidThermo<Thermo>::calce() const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName("e", this->group()),
            this->p_.mesh(),
            dimEnergy/dimMass
        )
    );
    volScalarField& psi(tPsi.ref());
    forAll(this->T_, celli)
    {
        const Thermo& mixture(volMixture_[celli]);
        psi[celli] = mixture.initializeEnergy
        (
            this->p_[celli],
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        );
    }
    forAll(this->T_.boundaryField(), patchi)
    {
        forAll(this->T_.boundaryField()[patchi], facei)
        {
            const Thermo& mixture(faceMixture_[patchi][facei]);
            psi.boundaryFieldRef()[patchi][facei] =
                mixture.initializeEnergy
                (
                    this->p_.boundaryField()[patchi][facei],
                    this->rho_.boundaryField()[patchi][facei],
                    this->e_.boundaryField()[patchi][facei],
                    this->T_.boundaryField()[patchi][facei]
                );
        }
    }
    return tPsi;
}

template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidThermo<Thermo>::calcT() const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName("T", this->group()),
            this->p_.mesh(),
            dimTemperature
        )
    );
    volScalarField& psi(tPsi.ref());
    forAll(this->T_, celli)
    {
        const Thermo& mixture(volMixture_[celli]);
        psi[celli] = mixture.TRhoE
        (
            this->T_[celli],
            this->rho_[celli],
            this->e_[celli]

        );
    }
    forAll(this->T_.boundaryField(), patchi)
    {
        forAll(this->T_.boundaryField()[patchi], facei)
        {
            const Thermo& mixture(faceMixture_[patchi][facei]);
            psi.boundaryFieldRef()[patchi][facei] =
                mixture.TRhoE
                (
                    this->T_.boundaryField()[patchi][facei],
                    this->rho_.boundaryField()[patchi][facei],
                    this->e_.boundaryField()[patchi][facei]
                );
        }
    }
    return tPsi;
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidThermo<Thermo>::TRhoE
(
    const scalarField& T,
    const scalarField& e,
    const label patchi
) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->T_.boundaryField()[patchi], facei)
    {
        const Thermo& mixture(faceMixture_[patchi][facei]);
        psi[facei] =
            mixture.TRhoE
            (
                T[facei],
                this->rho_.boundaryField()[patchi][facei],
                e[facei]
            );
    }
    return tPsi;
}


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidThermo<Thermo>::TRhoEi
(
    const scalar& T,
    const scalar& e,
    const label celli
) const
{
    return
        volMixture_[celli].TRhoE
        (
            T,
            this->rho_[celli],
            e
        );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidThermo<Thermo>::E() const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName("e", this->group()),
            this->p_.mesh(),
            dimEnergy/dimMass
        )
    );
    volScalarField& psi(tPsi.ref());
    forAll(this->T_, celli)
    {
        const Thermo& mixture(volMixture_[celli]);
        psi[celli] = mixture.Es
        (
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        );
    }
    forAll(this->T_.boundaryField(), patchi)
    {
        forAll(this->T_.boundaryField()[patchi], facei)
        {
            const Thermo& mixture(faceMixture_[patchi][facei]);
            psi.boundaryFieldRef()[patchi][facei] =
                mixture.Es
                (
                    this->rho_.boundaryField()[patchi][facei],
                    this->e_.boundaryField()[patchi][facei],
                    this->T_.boundaryField()[patchi][facei]
                );
        }
    }
    return tPsi;
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidThermo<Thermo>::e
(
    const volScalarField& rho,
    const volScalarField& e,
    const volScalarField& T
) const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName("e", this->group()),
            this->p_.mesh(),
            dimEnergy/dimMass
        )
    );
    volScalarField& psi(tPsi.ref());
    forAll(this->T_, celli)
    {
        const Thermo& mixture(volMixture_[celli]);
        psi[celli] = mixture.Es
        (
            rho[celli],
            e[celli],
            T[celli]
        );
    }
    forAll(this->T_.boundaryField(), patchi)
    {
        forAll(this->T_.boundaryField()[patchi], facei)
        {
            const Thermo& mixture(faceMixture_[patchi][facei]);
            psi.boundaryFieldRef()[patchi][facei] =
                mixture.Es
                (
                    rho.boundaryField()[patchi][facei],
                    e.boundaryField()[patchi][facei],
                    T.boundaryField()[patchi][facei]
                );
        }
    }
    return tPsi;
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidThermo<Thermo>::e
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->T_.boundaryField()[patchi], facei)
    {
        const Thermo& mixture(faceMixture_[patchi][facei]);
        psi[facei] =
            mixture.Es
            (
                rho[facei],
                e[facei],
                T[facei]
            );
    }
    return tPsi;
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidThermo<Thermo>::e
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const labelList& faceCells
) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(faceCells.size())
    );
    scalarField& psi = tPsi.ref();

    forAll(faceCells, facei)
    {
        const Thermo& mixture(volMixture_[faceCells[facei]]);
        psi[facei] =
            mixture.Es
            (
                rho[facei],
                e[facei],
                T[facei]
            );
    }
    return tPsi;
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidThermo<Thermo>::W() const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName("W", this->group()),
            this->p_.mesh(),
            dimMass/dimMoles
        )
    );
    volScalarField& psi(tPsi.ref());
    forAll(this->T_, celli)
    {
        const Thermo& mixture(volMixture_[celli]);
        psi[celli] = mixture.W();
    }
    forAll(this->T_.boundaryField(), patchi)
    {
        forAll(this->T_.boundaryField()[patchi], facei)
        {
            const Thermo& mixture(faceMixture_[patchi][facei]);
            psi.boundaryFieldRef()[patchi][facei] = mixture.W();
        }
    }
    return tPsi;
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidThermo<Thermo>::W(const label patchi) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->T_.boundaryField()[patchi], facei)
    {
        const Thermo& mixture(faceMixture_[patchi][facei]);
        psi[facei] = mixture.W();
    }
    return tPsi;
}


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidThermo<Thermo>::Wi(const label celli) const
{
    return volMixture_[celli].W();
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidThermo<Thermo>::Gamma() const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName("Gamma", this->group()),
            this->p_.mesh(),
            dimless
        )
    );
    volScalarField& psi(tPsi.ref());
    forAll(this->T_, celli)
    {
        const Thermo& mixture(volMixture_[celli]);
        psi[celli] = mixture.Gamma
        (
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        );
    }
    forAll(this->T_.boundaryField(), patchi)
    {
        forAll(this->T_.boundaryField()[patchi], facei)
        {
            const Thermo& mixture(faceMixture_[patchi][facei]);
            psi.boundaryFieldRef()[patchi][facei] =
                mixture.Gamma
                (
                    this->rho_.boundaryField()[patchi][facei],
                    this->e_.boundaryField()[patchi][facei],
                    this->T_.boundaryField()[patchi][facei]
                );
        }
    }
    return tPsi;
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidThermo<Thermo>::Gamma(const label patchi) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->T_.boundaryField()[patchi], facei)
    {
        const Thermo& mixture(faceMixture_[patchi][facei]);
        psi[facei] =
            mixture.Gamma
            (
                this->rho_.boundaryField()[patchi][facei],
                this->e_.boundaryField()[patchi][facei],
                this->T_.boundaryField()[patchi][facei]
            );
    }
    return tPsi;
}


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidThermo<Thermo>::Gammai(const label celli) const
{
    return volMixture_[celli].Gamma
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidThermo<Thermo>::Cp() const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName("Cp", this->group()),
            this->p_.mesh(),
            dimEnergy/dimMass/dimTemperature
        )
    );
    volScalarField& psi(tPsi.ref());
    forAll(this->T_, celli)
    {
        const Thermo& mixture(volMixture_[celli]);
        psi[celli] = mixture.Cp
        (
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        );
    }
    forAll(this->T_.boundaryField(), patchi)
    {
        forAll(this->T_.boundaryField()[patchi], facei)
        {
            const Thermo& mixture(faceMixture_[patchi][facei]);
            psi.boundaryFieldRef()[patchi][facei] =
                mixture.Cp
                (
                    this->rho_.boundaryField()[patchi][facei],
                    this->e_.boundaryField()[patchi][facei],
                    this->T_.boundaryField()[patchi][facei]
                );
        }
    }
    return tPsi;
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidThermo<Thermo>::Cp(const label patchi) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->T_.boundaryField()[patchi], facei)
    {
        const Thermo& mixture(faceMixture_[patchi][facei]);
        psi[facei] =
            mixture.Cp
            (
                this->rho_.boundaryField()[patchi][facei],
                this->e_.boundaryField()[patchi][facei],
                this->T_.boundaryField()[patchi][facei]
            );
    }
    return tPsi;
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidThermo<Thermo>::Cp
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->T_.boundaryField()[patchi], facei)
    {
        const Thermo& mixture(faceMixture_[patchi][facei]);
        psi[facei] =
            mixture.Cp
            (
                rho[facei],
                e[facei],
                T[facei]
            );
    }
    return tPsi;
}


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidThermo<Thermo>::Cpi(const label celli) const
{
    return volMixture_[celli].Cp
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidThermo<Thermo>::Cv() const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName("Cv", this->group()),
            this->p_.mesh(),
            dimEnergy/dimMass/dimTemperature
        )
    );
    volScalarField& psi(tPsi.ref());
    forAll(this->T_, celli)
    {
        const Thermo& mixture(volMixture_[celli]);
        psi[celli] = mixture.Cv
        (
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        );
    }
    forAll(this->T_.boundaryField(), patchi)
    {
        forAll(this->T_.boundaryField()[patchi], facei)
        {
            const Thermo& mixture(faceMixture_[patchi][facei]);
            psi.boundaryFieldRef()[patchi][facei] =
                mixture.Cv
                (
                    this->rho_.boundaryField()[patchi][facei],
                    this->e_.boundaryField()[patchi][facei],
                    this->T_.boundaryField()[patchi][facei]
                );
        }
    }
    return tPsi;
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidThermo<Thermo>::Cv(const label patchi) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->T_.boundaryField()[patchi], facei)
    {
        const Thermo& mixture(faceMixture_[patchi][facei]);
        psi[facei] =
            mixture.Cv
            (
                this->rho_.boundaryField()[patchi][facei],
                this->e_.boundaryField()[patchi][facei],
                this->T_.boundaryField()[patchi][facei]
            );
    }
    return tPsi;
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidThermo<Thermo>::Cv
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->T_.boundaryField()[patchi], facei)
    {
        const Thermo& mixture(faceMixture_[patchi][facei]);
        psi[facei] =
            mixture.Cv
            (
                rho[facei],
                e[facei],
                T[facei]
            );
    }
    return tPsi;
}


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidThermo<Thermo>::Cvi(const label celli) const
{
    return volMixture_[celli].Cv
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidThermo<Thermo>::CpByCv() const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName("CpByCv", this->group()),
            this->p_.mesh(),
            dimless
        )
    );
    volScalarField& psi(tPsi.ref());
    forAll(this->T_, celli)
    {
        const Thermo& mixture(volMixture_[celli]);
        psi[celli] = mixture.CpByCv
        (
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        );
    }
    forAll(this->T_.boundaryField(), patchi)
    {
        forAll(this->T_.boundaryField()[patchi], facei)
        {
            const Thermo& mixture(faceMixture_[patchi][facei]);
            psi.boundaryFieldRef()[patchi][facei] =
                mixture.CpByCv
                (
                    this->rho_.boundaryField()[patchi][facei],
                    this->e_.boundaryField()[patchi][facei],
                    this->T_.boundaryField()[patchi][facei]
                );
        }
    }
    return tPsi;
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidThermo<Thermo>::CpByCv(const label patchi) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->T_.boundaryField()[patchi], facei)
    {
        const Thermo& mixture(faceMixture_[patchi][facei]);
        psi[facei] =
            mixture.CpByCv
            (
                this->rho_.boundaryField()[patchi][facei],
                this->e_.boundaryField()[patchi][facei],
                this->T_.boundaryField()[patchi][facei]
            );
    }
    return tPsi;
}


template<class Thermo>
Foam::tmp<Foam::scalarField>
Foam::multicomponentFluidThermo<Thermo>::CpByCv
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();

    forAll(this->T_.boundaryField()[patchi], facei)
    {
        const Thermo& mixture(faceMixture_[patchi][facei]);
        psi[facei] =
            mixture.CpByCv
            (
                rho[facei],
                e[facei],
                T[facei]
            );
    }
    return tPsi;
}

// ************************************************************************* //
