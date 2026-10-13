/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "reactingDiameterModel.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(reactingDiameterModel, 0);
    addToRunTimeSelectionTable(diameterModel, reactingDiameterModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::reactingDiameterModel::reactingDiameterModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    diameter(mesh, dict, phaseName),
    rate_(diameterReactionRate::New(*this, dict)),
    pName_(dict.lookupOrDefault("pName", word("p"))),
    TName_(dict.lookupOrDefault("TName", IOobject::groupName("T", phaseName))),
    oxidiserPhase_(dict.lookupOrDefault("oxidiserPhase", word::null)),
    oxidiserSpecies_
    (
        dict.lookupOrDefault("oxidiserSpecies", wordList())
    ),
    YoxidiserRef_
    (
        oxidiserSpecies_.size()
      ? dict.lookup<scalar>("YoxidiserRef")
      : 1.0
    ),
    N_
    (
        IOobject
        (
            IOobject::groupName("N", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", inv(dimVolume), 1.0)
    ),
    d0_
    (
        IOobject
        (
            IOobject::groupName("d0", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->d_
    ),
    dVdt_
    (
        IOobject
        (
            IOobject::groupName("dVdt", phaseName),
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimVolume/dimTime, 0.0)
    )
{
    this->requireD();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::reactingDiameterModel::~reactingDiameterModel()
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::reactingDiameterModel::oxidiserFactor() const
{
    const fvMesh& mesh = this->d_.mesh();
    const volScalarField& alpha =
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", this->d_.group())
        );
    const volScalarField& alphaOx =
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", oxidiserPhase_)
        );

    // Fraction of the fluid around the particles that carries the oxidiser
    tmp<volScalarField> tf
    (
        volScalarField::New
        (
            IOobject::groupName("oxidiserFactor", this->d_.group()),
            min
            (
                max(alphaOx, scalar(0))
               /max(1.0 - alpha, dimensionedScalar(dimless, small)),
                scalar(1)
            )
        )
    );

    if (oxidiserSpecies_.size())
    {
        volScalarField Yox
        (
            volScalarField::New
            (
                "Yoxidiser",
                mesh,
                dimensionedScalar(dimless, 0)
            )
        );
        forAll(oxidiserSpecies_, i)
        {
            Yox +=
                mesh.lookupObject<volScalarField>
                (
                    IOobject::groupName(oxidiserSpecies_[i], oxidiserPhase_)
                );
        }
        tf.ref() *= min(max(Yox, scalar(0))/YoxidiserRef_, scalar(1));
    }

    return tf;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::reactingDiameterModel::update()
{
    const volScalarField& alpha =
        this->d_.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", this->d_.group())
        );

    forAll(this->d_, celli)
    {
        if (alpha[celli] <= 1e-4)
        {
            this->d_[celli] = d0_[celli];
        }
        else
        {
            this->d_[celli] = Foam::max(1e-10, Foam::min(this->d_[celli], d0_[celli]));
        }
    }
    this->d_.correctBoundaryConditions();

    const dimensionedScalar piBy6
    (
        "piBy6",
        dimless,
        constant::mathematical::pi/6.0
    );
    N_ =
        alpha
       /max(piBy6*pow3(d0_), dimensionedScalar(dimVolume, small));
}


void Foam::diameterModels::reactingDiameterModel::solve()
{
    const volScalarField& p
    (
        this->d_.mesh().lookupObject<volScalarField>(pName_)
    );
    const volScalarField& T
    (
        this->d_.mesh().lookupObject<volScalarField>(TName_)
    );
    solve(p, T);
}


void Foam::diameterModels::reactingDiameterModel::solve
(
    const volScalarField& pi,
    const volScalarField& T
)
{
    const volScalarField& alpha =
        this->d_.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", this->d_.group())
        );

    volScalarField VOld
    (
        IOobject::groupName("V", this->d_.group()),
        this->V()
    );
    this->storeAndBlendOld(VOld);

    const volScalarField dOld(this->d_);
    this->storeAndBlendOld(this->d_);

    volScalarField dDdt(-2.0*rate_->dDdt(pi, T));
    if (!oxidiserPhase_.empty())
    {
        dDdt *= oxidiserFactor();
    }
    this->blendDelta(dDdt);

    const dimensionedScalar& dT(this->d_.time().deltaT());

    //- Change dRdt to dDdt
    this->d_ += dT*dDdt;
    forAll(this->d_, celli)
    {
        if (alpha[celli] <= 1e-4)
        {
            this->d_[celli] = d0_[celli];
        }
        else
        {
            this->d_[celli] = Foam::max(1e-10, Foam::min(this->d_[celli], d0_[celli]));
        }
    }
    this->d_.correctBoundaryConditions();

    dVdt_ = (this->V() - VOld)/dT;
    dDdt = (this->d_ - dOld)/dT;

    //- Compute actual delta for the time step knowing the blended
    dDdt = this->calcAndStoreDelta(dDdt);
    dVdt_ = this->calcAndStoreDelta(dVdt_);

    const surfaceScalarField& phi
    (
        this->d_.mesh().lookupObject<surfaceScalarField>
        (
            IOobject::groupName("phi", this->d_.group())
        )
    );

    volScalarField deltaPhid(fvc::div(phi, d_) - d_*fvc::div(phi));
    this->storeAndBlendDelta(deltaPhid);
    this->d_ = dOld - dT*(deltaPhid - dDdt);
    forAll(this->d_, celli)
    {
        if (alpha[celli] <= 1e-4)
        {
            this->d_[celli] = d0_[celli];
        }
        else
        {
            this->d_[celli] = Foam::max(1e-10, Foam::min(this->d_[celli], dOld[celli]));
        }
    }
    this->d_.correctBoundaryConditions();

    dVdt_ = (this->V() - VOld)/dT;
    dVdt_.min(0.0);
    dDdt = (this->d_ - dOld)/dT;
    dDdt.min(0.0);

    dDdt = this->calcAndStoreDelta(dDdt);
    dVdt_ = this->calcAndStoreDelta(dVdt_);
}


Foam::tmp<Foam::volScalarField>
Foam::diameterModels::reactingDiameterModel::dMdt() const
{
    const volScalarField& rho
    (
        this->d_.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("rho", this->d_.group())
        )
    );

    return dVdt_*rho*N_;
}

// ************************************************************************* //
