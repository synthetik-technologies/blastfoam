/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "compressibilityCorrection.H"
#include "fluidThermo.H"
#include "fluidBlastThermo.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressible::correction::correction
(
    dictionary& dict
)
:
    XiStar_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "XiStar",
            dict,
            2.0
        )
    ),

    Mat0_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Mat0",
            dict,
            0.25
        )
    ),

    rProdLim_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "rProdLim",
            dict,
            20.0
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::compressible::correction::read(const dictionary& dict)
{
    XiStar_.readIfPresent(dict);
    Mat0_.readIfPresent(dict);
    return true;
}


Foam::tmp<Foam::volScalarField::Internal> Foam::compressible::correction::Mt
(
    const volScalarField& k
)
{
    const basicThermo& thermo = k.mesh().lookupObject<basicThermo>
    (
        IOobject::groupName
        (
            physicalProperties::typeName,
            k.group()
        )
    );
    if (isA<fluidBlastThermo>(thermo))
    {
        return
            sqrt(2.0*k())
           /(dynamicCast<const fluidBlastThermo>(thermo).speedOfSound()());
    }
    else if (isA<fluidThermo>(thermo))
    {
        const fluidThermo& fluid = dynamicCast<const fluidThermo>(thermo);
        return sqrt(2.0*k()/(fluid.Cp()()/fluid.Cv()()/fluid.psi()()));
    }
    else
    {
        FatalErrorInFunction
            << "Only fluidThermo thermodynamic models can be used with "
            << "compressibility corrections" << endl
            << abort(FatalError);
    }
    return tmp<volScalarField::Internal>();
}


Foam::tmp<Foam::volScalarField::Internal> Foam::compressible::correction::MtSqr
(
    const volScalarField& k
)
{
    const basicThermo& thermo = k.mesh().lookupObject<basicThermo>
    (
        IOobject::groupName
        (
            physicalProperties::typeName,
            k.group()
        )
    );
    if (isA<fluidBlastThermo>(thermo))
    {
        return
            2.0*k()
           /sqr(dynamicCast<const fluidBlastThermo>(thermo).speedOfSound()());
    }
    else if (isA<fluidThermo>(thermo))
    {
        const fluidThermo& fluid = dynamicCast<const fluidThermo>(thermo);
        return 2.0*k()/(fluid.Cp()()/fluid.Cv()()/fluid.psi()());
    }
    else
    {
        FatalErrorInFunction
            << "Only fluidThermo thermodynamic models can be used with "
            << "compressibility corrections" << endl
            << abort(FatalError);
    }
    return tmp<volScalarField::Internal>();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::compressible::correction::MtSqrByk
(
    const volScalarField& k
)
{
    const basicThermo& thermo = k.mesh().lookupObject<basicThermo>
    (
        IOobject::groupName
        (
            physicalProperties::typeName,
            k.group()
        )
    );
    if (isA<fluidBlastThermo>(thermo))
    {
        return
            2.0
            /sqr(dynamicCast<const fluidBlastThermo>(thermo).speedOfSound()());
    }
    else if (isA<fluidThermo>(thermo))
    {
        const fluidThermo& fluid = dynamicCast<const fluidThermo>(thermo);
        return 2.0/(fluid.Cp()()/fluid.Cv()()/fluid.psi()());
    }
    else
    {
        FatalErrorInFunction
            << "Only fluidThermo thermodynamic models can be used with "
            << "compressibility corrections" << endl
            << abort(FatalError);
    }
    return tmp<volScalarField::Internal>();
}

Foam::tmp<Foam::volScalarField::Internal> Foam::compressible::correction::Fcorr
(
    const volScalarField& k
) const
{
    volScalarField::Internal Mat(Mt(k));
    return (sqr(Mat) - sqr(Mat0_))*pos(Mat - Mat0_);

    // alphaSarkar_ = 0.5;
    // return (1.0/(1.0 + alphaSarkar_*sqr(Mat));

}

void Foam::compressible::correction::correct
(
    const volScalarField& k,
    const dimensionedScalar& beta0,
    const dimensionedScalar& betaStar0,
    tmp<volScalarField::Internal>& beta,
    tmp<volScalarField::Internal>& betaStar
) const
{
    const volScalarField::Internal f(Fcorr(k));
    betaStar = betaStar0*(1.0 + XiStar_*f);
    beta = beta0 - betaStar()*XiStar_*Fcorr(k);
}


void Foam::compressible::correction::correct
(
    const volScalarField& k,
    const dimensionedScalar& betaStar0,
    tmp<volScalarField::Internal>& beta,
    tmp<volScalarField::Internal>& betaStar
) const
{
    const volScalarField::Internal f(Fcorr(k));
    betaStar = betaStar0*(1.0 + XiStar_*f);
    beta = beta - betaStar()*XiStar_*Fcorr(k);
}


void Foam::compressible::correction::limitG
(
    volScalarField::Internal& G,
    const volScalarField::Internal& k,
    const volScalarField::Internal& omega,
    const volScalarField::Internal& betaStar
) const
{
    G = min(G, rProdLim_*betaStar*k*omega);
}


// ************************************************************************* //
