/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "errorType.H"
#include "IFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(errorType, 0);
    defineRunTimeSelectionTable(errorType, dictionary);

    template<>
    const char* NamedEnum<errorType::Reduction, 7>::names[] =
    {
        "min",
        "minMagSqr",
        "max",
        "maxMagSqr",
        "average",
        "weightedAverage",
        "latest"
    };
    const NamedEnum<errorType::Reduction, 7>
        errorType::timeReductionTypeNames;


    template<>
    const char* NamedEnum<errorType::Target, 4>::names[] =
    {
        "zero",
        "target",
        "lowerBound",
        "upperBound"
    };

    const NamedEnum<errorType::Target, 4> errorType::targetTypeNames;
}


Foam::label Foam::errorType::getConfigurationNo()
{
    string str(getEnv("BLAST_CONFIG_NO"));
    if (str.size())
    {
        IStringStream iss(str);
        return readLabel(iss);
    }
    return 0;
}

Foam::label Foam::errorType::getNConfigurations()
{
    string str(getEnv("N_BLAST_CONFIGS"));
    if (str.size())
    {
        IStringStream iss(str);
        return readLabel(iss);
    }
    return 0;
}

const Foam::label Foam::errorType::configuration =
    ::Foam::errorType::getConfigurationNo();

const Foam::label Foam::errorType::nConfigurations =
    ::Foam::errorType::getNConfigurations();


Foam::scalar Foam::errorType::reduceValue(const scalar f)
{
    switch (timeReduction_)
    {
        case MIN:
        {
            return min(value_, f);
        }
        case MIN_MAG_SQR:
        {
            return magSqr(value_) < magSqr(f) ? value_ : f;
        }
        case MAX:
        {
            return max(value_, f);
        }
        case MAX_MAG_SQR:
        {
            return magSqr(value_) > magSqr(f) ? value_ : f;
        }
        case AVERAGE:
        {
            intValue_ += f;
            intWeight_ += 1.0;
            return intValue_/intWeight_;
        }
        case WEIGHTED_AVERAGE:
        {
            intValue_ += f*runTime_.deltaTValue();
            intWeight_ += runTime_.deltaTValue();
            return intValue_/intWeight_;
        }
        case LATEST:
        {
            return f;
        }
    }
    return intValue_/intWeight_;
}


Foam::scalar Foam::errorType::reduceField
(
    const label method,
    const scalarField& f,
    const scalarField& w
)
{
    switch (method)
    {
        case MIN:
        {
            return gMin(f);
        }
        case MIN_MAG_SQR:
        {
            return gMinMagSqr(f);
        }
        case MAX:
        {
            return gMax(f);
        }
        case MAX_MAG_SQR:
        {
            return gMaxMagSqr(f);
        }
        case AVERAGE:
        {
            return gAverage(f);
        }
        case WEIGHTED_AVERAGE:
        {
            return gSum(f*w)/gSum(w);
        }
        default:
        {
            FatalErrorInFunction
                << "Unsupported field reduction method" << endl
                << abort(FatalError);
        }
    }
    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorType::errorType
(
    const Time& runTime,
    const dictionary& dict,
    const word& region
)
:
    errorType
    (
        runTime,
        dict,
        region,
        timeReductionTypeNames.read(dict.lookup("timeReduction"))
    )
{}

Foam::errorType::errorType
(
    const Time& runTime,
    const dictionary& dict,
    const word& region,
    const Reduction& timeReduction
)
:
    name_(dict.dictName()),
    runTime_(runTime),
    regionName_(region),
    normalise_(dict.lookupOrDefault("normalise", true)),
    value_(0.0),

    timeReduction_(timeReduction),
    intValue_(0.0),
    intWeight_(0.0),

    targetType_(targetTypeNames.read(dict.lookup("targetType"))),
    targetValue_
    (
        targetType_ != ZERO
      ? readConfigValue<scalar>("targetValue", dict)
      : 0.0
    ),
    weight_(dict.lookupOrDefault<scalar>("weight", 1.0))
{
    clear();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorType::~errorType()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::scalar Foam::errorType::error() const
{
    if
    (
        (targetType_ == LOWER_BOUND && value_ > targetValue_)
     || (targetType_ == UPPER_BOUND && value_ < targetValue_)
    )
    {
        return 0.0;

    }
    if (mag(targetValue_) > small)
    {
        if (normalise_)
        {
            return mag(value_ - targetValue_)/targetValue_;
        }
        else
        {
            return mag(value_ - targetValue_);
        }
    }
    return mag(value_);
}


void Foam::errorType::clear()
{
    switch (timeReduction_)
    {
        case MIN:
        case MIN_MAG_SQR:
        {
            value_ = great;
            return;
        }
        case MAX:
        case MAX_MAG_SQR:
        {
            value_ = -great;
            return;
        }
        default:
        {
            value_ = 0.0;
            return;
        }
    }
}


// ************************************************************************* //
