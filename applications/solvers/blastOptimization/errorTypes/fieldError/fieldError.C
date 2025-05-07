/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "fieldError.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace errorTypes
{
    defineTypeNameAndDebug(field, 0);
    addToRunTimeSelectionTable(errorType, field, dictionary);

    template<>
    const char* NamedEnum<field::Operation, 5>::names[] =
    {
        "none",
        "laplacian",
        "ddt",
        "magGrad",
        "function"
    };
    const NamedEnum<field::Operation, 5> field::operationNames_;

    template<>
    const char* NamedEnum<field::FieldReduction, 7>::names[] =
    {
        "min",
        "minMagSqr",
        "max",
        "maxMagSqr",
        "average",
        "weightedAverage",
        "probe"
    };

    const NamedEnum<field::FieldReduction, 7>
        field::fieldReductionTypeNames;
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorTypes::field::field
(
    const Time& runTime,
    const dictionary& dict,
    const word& region
)
:
    field(runTime, dict, region, operationNames_.read(dict.lookup("operator")))
{}


Foam::errorTypes::field::field
(
    const Time& runTime,
    const dictionary& dict,
    const word& region,
    const Operation op
)
:
    errorType(runTime, dict, region),
    fieldName_(dict.lookup("field")),
    operation_(op),
    fieldReduction_
    (
        fieldReductionTypeNames.read(dict.lookup("fieldReduction"))
    ),
    timeVarying_(false)
{
    if (operation_ == FUNCTION)
    {
        const polyMesh& mesh = this->mesh();
        timeVarying_ = dict.lookup<Switch>("timeVarying");
        List<word> cmptNames;
        forAll(mesh.geometricD(), i)
        {
            if (mesh.geometricD()[i] > 0)
            {
                cmpts_.append(i);
                cmptNames.append(vector::componentNames[i]);
            }
        }
        if (timeVarying_)
        {
            cmpts_.append(-1);
        }
        if (cmpts_.size() == 1)
        {
            Info<< "Reading Function1<scalar>" << endl;
            func1_ = Function1<scalar>::New
            (
                "function",
                runTime.userUnits(),
                dimless,
                dict
            );
        }
        else if (cmpts_.size() == 2)
        {
            Info<< "Reading Function2<scalar>" << endl;
            func2_ = Function2<scalar>::New
            (
                "function",
                dimless,
                dimless,
                dimless,
                dict
            );
        }
        else if (cmpts_.size() == 3)
        {
            Info<< "Reading Function2<scalar>" << endl;
            func3_ = Function3<scalar>::New
            (
                "function",
                dimless,
                dimless,
                dimless,
                dimless,
                dict
            );
        }
        // else if (cmpts_.size() == 4)
        // {
        //     func4_ = Function4<scalar>::New("function", dict);
        // }
    }

    if (fieldReduction_ == PROBE)
    {
        probeLocation_ = dict.lookup<vector>("location");
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorTypes::field::~field()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::errorTypes::field::performOperation
(
    const Foam::volScalarField& f
) const
{
    switch (operation_)
    {
        case NONE:
        {
            return f;
        }
        case DDT:
        {
            return fvc::ddt(f);
        }
        case LAPLACIAN:
        {
            return fvc::laplacian(f);
        }
        case MAG_GRAD:
        {
            return mag(fvc::grad(f));
        }
        case FUNCTION:
        {
            tmp<volScalarField> ttarget =
                volScalarField::New
                (
                    "target",
                    mesh(),
                    dimensionedScalar(f.dimensions(), 0.0)
                );
            volScalarField& target = ttarget.ref();
            const volVectorField& C = mesh().C();
            const scalar& t = mesh().time().value();
            if (func1_.valid())
            {
                forAll(C, i)
                {
                    target[i] = func1_->value
                    (
                        cmpts_[0] < 0 ? t : C[i][cmpts_[0]]
                    );
                }
            }
            else if (func2_.valid())
            {
                forAll(C, i)
                {
                    target[i] = func2_->value
                    (
                        C[i][cmpts_[0]],
                        cmpts_[1] < 0 ? t : C[i][cmpts_[1]]
                    );
                }
            }
            else if (func3_.valid())
            {
                forAll(C, i)
                {
                    target[i] = func3_->value
                    (
                        C[i][cmpts_[0]],
                        C[i][cmpts_[1]],
                        cmpts_[2] < 0 ? t : C[i][cmpts_[2]]
                    );
                }
            }
            else
            {
                forAll(C, i)
                {
                    target[i] = func4_->value
                    (
                        C[i][cmpts_[0]],
                        C[i][cmpts_[1]],
                        C[i][cmpts_[2]],
                        t
                    );
                }
            }
            return target - f;
        }
    }
    return tmp<volScalarField>();
}


void Foam::errorTypes::field::update()
{
    // Look up field
    const volScalarField& field
    (
        this->mesh().lookupObject<volScalarField>(fieldName_)
    );

    // Perform operation to field
    tmp<volScalarField> topField(performOperation(field));

    if (fieldReduction_ == PROBE)
    {
        const label celli = mesh().findCell(probeLocation_);
        scalar val = -great;
        if (celli >= 0)
        {
            val = topField()[celli];
        }
        reduce(val, maxOp<scalar>());
        value_ = reduceValue(val);
    }
    else
    {
        // Reduce the field and update time based values
        value_ = reduceValue
        (
            reduceField(fieldReduction_, topField(), this->mesh().V())
        );
    }
}

// ************************************************************************* //
