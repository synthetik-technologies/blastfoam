/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "SpalartAllmaras_comp.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "fluxScheme.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField>
SpalartAllmaras_comp<BasicMomentumTransportModel>::chi() const
{
    return volScalarField::New(typedName("chi"), nuTilda_/this->nu());
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
SpalartAllmaras_comp<BasicMomentumTransportModel>::fv1
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(typedName("chi3"), pow3(chi));
    return volScalarField::New(typedName("fv1"), chi3/(chi3 + pow3(Cv1_)));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
SpalartAllmaras_comp<BasicMomentumTransportModel>::fv2
(
    const volScalarField::Internal& chi,
    const volScalarField::Internal& fv1
) const
{
    return volScalarField::Internal::New
    (
        typedName("fv2"),
        1.0 - chi/(1.0 + chi*fv1)
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
SpalartAllmaras_comp<BasicMomentumTransportModel>::ft2
(
    const volScalarField::Internal& chi
) const
{
    return volScalarField::Internal::New
    (
        typedName("ft2"),
        Ct3_*exp(-Ct4_*sqr(chi))
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
SpalartAllmaras_comp<BasicMomentumTransportModel>::Stilda
(
    const volTensorField::Internal& gradU,
    const volScalarField::Internal& chi,
    const volScalarField::Internal& fv1
) const
{
    const volScalarField::Internal Omega
    (
        typedName("Omega"),
        ::sqrt(2.0)*mag(skew(gradU))
    );

    return volScalarField::Internal::New
    (
        typedName("Stilda"),
        (
            max
            (
                Omega
              + fv2(chi, fv1)*nuTilda_()/sqr(kappa_*this->y()()),
                Cs_*Omega
            )
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
SpalartAllmaras_comp<BasicMomentumTransportModel>::fw
(
    const volScalarField::Internal& Stilda
) const
{
    const volScalarField::Internal r
    (
        typedName("r"),
        min
        (
            nuTilda_()
           /(
               max
               (
                   Stilda,
                   dimensionedScalar(Stilda.dimensions(), small)
               )
              *sqr(kappa_*this->y()())
            ),
            scalar(10.0)
        )
    );

    const volScalarField::Internal g
    (
        typedName("g"),
        r + Cw2_*(pow6(r) - r)
    );

    return volScalarField::Internal::New
    (
        typedName("fw"),
        g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0)
    );
}


template<class BasicMomentumTransportModel>
void SpalartAllmaras_comp<BasicMomentumTransportModel>::correctNut
(
    const volScalarField& fv1
)
{
    this->nut_ = nuTilda_*fv1;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
void SpalartAllmaras_comp<BasicMomentumTransportModel>::correctNut()
{
    correctNut(fv1(this->chi()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
SpalartAllmaras_comp<BasicMomentumTransportModel>::SpalartAllmaras_comp
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),
    compressible::correction
    (
        this->coeffDict_,
        MODEL::OTHER,
        *this
    ),

    sigmaNut_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaNut",
            this->coeffDict_,
            2.0/3.0
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),

    Cb1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb1",
            this->coeffDict_,
            0.1355
        )
    ),
    Cb2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb2",
            this->coeffDict_,
            0.622
        )
    ),
    Cw1_(Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_),
    Cw2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            this->coeffDict_,
            0.3
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw3",
            this->coeffDict_,
            2.0
        )
    ),
    Cv1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cv1",
            this->coeffDict_,
            7.1
        )
    ),
    Ct3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ct3",
            this->coeffDict_,
            1.2
        )
    ),
    Ct4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ct4",
            this->coeffDict_,
            0.5
        )
    ),
    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.3
        )
    ),
    C5_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C5",
            this->coeffDict_,
            0.0 //3.5
        )
    ),

    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    if (type == typeName || type == baseName())
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool SpalartAllmaras_comp<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        sigmaNut_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());

        Cb1_.readIfPresent(this->coeffDict());
        Cb2_.readIfPresent(this->coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_;
        Cw2_.readIfPresent(this->coeffDict());
        Cw3_.readIfPresent(this->coeffDict());
        Cv1_.readIfPresent(this->coeffDict());
        Ct3_.readIfPresent(this->coeffDict());
        Ct4_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());
        C5_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
SpalartAllmaras_comp<BasicMomentumTransportModel>::DnuTildaEff() const
{
    return volScalarField::New
    (
        "DnuTildaEff",
        (nuTilda_ + this->nu())/sigmaNut_
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> SpalartAllmaras_comp<BasicMomentumTransportModel>::k() const
{
    return volScalarField::New
    (
        "k",
        this->mesh_,
        dimensionedScalar(dimensionSet(0, 2, -2, 0, 0), 0)
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
SpalartAllmaras_comp<BasicMomentumTransportModel>::epsilon() const
{
    WarningInFunction
        << "Turbulence kinetic energy dissipation rate not defined for "
        << "Spalart-Allmaras model. Returning zero field"
        << endl;

    return volScalarField::New
    (
        "epsilon",
        this->mesh_,
        dimensionedScalar(dimensionSet(0, 2, -3, 0, 0), 0)
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
SpalartAllmaras_comp<BasicMomentumTransportModel>::omega() const
{
    WarningInFunction
        << "Turbulence specific dissipation rate not defined for "
        << "Spalart-Allmaras model. Returning zero field"
        << endl;

    return volScalarField::New
    (
        "omega",
        this->mesh_,
        dimensionedScalar(dimless/dimTime, 0)
    );
}


template<class BasicMomentumTransportModel>
void SpalartAllmaras_comp<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

    const volVectorField gradNuTilda(fvc::grad(nuTilda_));
    const volTensorField gradU(fvc::grad(this->U()));
    const volScalarField::Internal Stilda(this->Stilda(gradU, chi, fv1));
    const volScalarField::Internal ft2(this->ft2(chi));

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(alpha, rho, nuTilda_)
      + fvm::div(alphaRhoPhi, nuTilda_)
      - fvm::laplacian(alpha*rho*DnuTildaEff(), nuTilda_)
      - alpha*rho*(Cb2_/sigmaNut_)*magSqr(gradNuTilda)
     ==
        alpha()*rho()*Cb1_*(1.0 - ft2)*Stilda*nuTilda_()
      - fvm::Sp
        (
            alpha()*rho()
           *(Cw1_*fw(Stilda) - (Cb1_/sqr(kappa_))*ft2)
           *nuTilda_()/sqr(this->y()()),
            nuTilda_
        )

        // Compressiblity correction
      - (this->nu() + nuTilda_)/sigmaNut_
       *(fvc::grad(rho) & gradNuTilda)

      + fvModels.source(alpha, rho, nuTilda_)
    );

    if (C5_.value() > small)
    {
        nuTildaEqn.ref() +=
            fvm::Sp
            (
                C5_*alpha()*rho()*nuTilda_()*(gradU() && gradU())
               /sqr(this->speedOfSound()),
                nuTilda_
            );
    }

    nuTildaEqn.ref().relax();
    fvConstraints.constrain(nuTildaEqn.ref());
    Foam::solve(nuTildaEqn);
    fvConstraints.constrain(nuTilda_);
    bound(nuTilda_, dimensionedScalar(nuTilda_.dimensions(), 0));
    nuTilda_.correctBoundaryConditions();

    correctNut(fv1);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
