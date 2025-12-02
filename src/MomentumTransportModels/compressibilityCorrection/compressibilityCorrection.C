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


template<>
const char* Foam::NamedEnum
<
    Foam::compressible::correction::F_CC_Method,
    Foam::compressible::correction::F_CC_Method::SIZE_F
>::names[] =
{
    "none",
    "Wilcox",
    "Zeman"
};
const Foam::NamedEnum
<
    Foam::compressible::correction::F_CC_Method,
    Foam::compressible::correction::F_CC_Method::SIZE_F
> Foam::compressible::correction::F_CC_MethodNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressible::correction::correction
(
    dictionary& dict,
    const MODEL model,
    const momentumTransportModel& turb,
    const F_CC_Method ccMethod
)
:
    turb_(turb),

    FMethod_(ccMethod),
    XiStar_("XiStar", dimless, 2.0),
    lambda_("lambda", dimless, 0.66),
    Mt0_("Mt0", dimless, 0.25),

    pressureDilation_(false),
    alphaK2_("alphaK2", dimless, 0.15),
    alphaK3_("alphaK3", dimless, 0.2)
{
    if (model == MODEL::K_OMEGA)
    {
        FMethod_ = F_CC_MethodNames_
        [
            dict.lookupOrAddDefault<word>
            (
                "compressibilityCorrection",
                F_CC_MethodNames_[FMethod_]
            )
        ];

        switch (FMethod_)
        {
            case F_CC_Method::NONE:
            {
                break;
            }
            case F_CC_Method::WILCOX:
            {
                XiStar_ =
                    dimensioned<scalar>::lookupOrAddToDict
                    (
                        "XiStar",
                        dict,
                        2.0
                    );
                Mt0_ =
                    dimensioned<scalar>::lookupOrAddToDict
                    (
                        "Mt0",
                        dict,
                        0.25
                    );
                break;
            }
            case F_CC_Method::ZEMAN:
            {
                XiStar_ =
                    dimensioned<scalar>::lookupOrAddToDict
                    (
                        "XiStar",
                        dict,
                        0.75
                    );
                lambda_ =
                    dimensioned<scalar>::lookupOrAddToDict
                    (
                        "lambda",
                        dict,
                        0.66
                    );
                Mt0_ =
                    dimensioned<scalar>::lookupOrAddToDict
                    (
                        "Mt0",
                        dict,
                        0.2
                    );
                break;
            }
            default:
            {
                NotImplemented;
            }
        }
    }
    else if (model == MODEL::K_EPSILON)
    {
        pressureDilation_ =
            dict.lookupOrAddDefault("pressureDilation", false);
    }


    if (pressureDilation_)
    {
        alphaK2_ =
            dimensioned<scalar>::lookupOrAddToDict
            (
                "alphaK2",
                dict,
                0.15
            );
        alphaK3_ =
            dimensioned<scalar>::lookupOrAddToDict
            (
                "alphaK3",
                dict,
                0.3
            );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::compressible::correction::read(const dictionary& dict)
{
    if (dict.found("compressibilityCorrection"))
    {
        FMethod_ =
            F_CC_MethodNames_
            [
                dict.lookup("compressibilityCorrection")
            ];
    }
    XiStar_.readIfPresent(dict);
    lambda_.readIfPresent(dict);
    Mt0_.readIfPresent(dict);

    dict.readIfPresent
    (
        "pressureDilation",
        pressureDilation_
    );

    alphaK2_.readIfPresent(dict);
    alphaK3_.readIfPresent(dict);

    return true;

}

Foam::tmp<Foam::volScalarField::Internal>
Foam::compressible::correction::speedOfSound() const
{
    const viscosity& visc = turb_.properties();
    if (isA<fluidBlastThermo>(visc))
    {
        return dynamicCast<const fluidBlastThermo>(visc).speedOfSound()();
    }
    else if (isA<fluidThermo>(visc))
    {
        const fluidThermo& fluid = dynamicCast<const fluidThermo>(visc);
        return sqrt(fluid.Cp()()/fluid.Cv()()/fluid.psi()());
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
Foam::compressible::correction::Mt() const
{
    const viscosity& visc = turb_.properties();
    tmp<volScalarField> tk = turb_.k();
    const volScalarField& k = tk();
    if (isA<fluidBlastThermo>(visc))
    {
        return
            sqrt(2.0*k())
           /(dynamicCast<const fluidBlastThermo>(visc).speedOfSound()());
    }
    else if (isA<fluidThermo>(visc))
    {
        const fluidThermo& fluid = dynamicCast<const fluidThermo>(visc);
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


Foam::tmp<Foam::volScalarField::Internal> Foam::compressible::correction::MtSqr() const
{
    const viscosity& visc = turb_.properties();
    tmp<volScalarField> tk = turb_.k();
    const volScalarField& k = tk();
    if (isA<fluidBlastThermo>(visc))
    {
        return
            2.0*k()
           /sqr(dynamicCast<const fluidBlastThermo>(visc).speedOfSound()());
    }
    else if (isA<fluidThermo>(visc))
    {
        const fluidThermo& fluid = dynamicCast<const fluidThermo>(visc);
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
Foam::compressible::correction::MtSqrByk() const
{
    const viscosity& visc = turb_.properties();
    if (isA<fluidBlastThermo>(visc))
    {
        return
            2.0
            /sqr(dynamicCast<const fluidBlastThermo>(visc).speedOfSound()());
    }
    else if (isA<fluidThermo>(visc))
    {
        const fluidThermo& fluid = dynamicCast<const fluidThermo>(visc);
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


Foam::tmp<Foam::volScalarField::Internal>
Foam::compressible::correction::Fcorr() const
{
    tmp<volScalarField::Internal> tMt;
    return Fcorr(tMt);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::compressible::correction::Fcorr
(
    tmp<volScalarField::Internal>& tMt
) const
{
    if (FMethod_ == F_CC_Method::NONE)
    {
        return volScalarField::Internal::New
        (
            "Fcorr",
            turb_.mesh(),
            0.0
        );
    }
    if (!tMt.valid())
    {
        tMt = this->Mt();
    }
    return Fcorr(tMt());
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::compressible::correction::Fcorr
(
    const volScalarField::Internal& Mt
) const
{
    switch (FMethod_)
    {
        case F_CC_Method::NONE:
        {
            return volScalarField::Internal::New
            (
                "Fcorr",
                Mt.mesh(),
                0.0
            );
        }
        case F_CC_Method::WILCOX:
        {
            return max(sqr(Mt) - sqr(Mt0_), 0.0);
        }
        case F_CC_Method::ZEMAN:
        {
            return max(1.0 - exp(-sqr((Mt - Mt0_)/lambda_)), 0.0);
        }
        default:
        {
            NotImplemented;
            return tmp<volScalarField::Internal>();
        }
    }
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::compressible::correction::pressureDialationSource
(
    const volScalarField::Internal& G
) const
{
    tmp<volScalarField::Internal> tMt;
    return pressureDialationSource(G, tMt);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::compressible::correction::pressureDialationSource
(
    const volScalarField::Internal& G,
    tmp<volScalarField::Internal>& tMt
) const
{
    if (!pressureDilation_)
    {
        return volScalarField::Internal::New
        (
            "pressureDilation",
            G.mesh(),
            dimensionedScalar(G.dimensions(), 0.0)
        );
    }

    if (!tMt.valid()) tMt = this->Mt();
    return sqr(tMt())*(-alphaK2_*G + alphaK3_*turb_.epsilon()()());
}


// ************************************************************************* //
