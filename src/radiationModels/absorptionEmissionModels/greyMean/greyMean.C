/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "greyMean.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(greyMean, 0);

    addToRunTimeSelectionTable
    (
        absorptionEmissionModel,
        greyMean,
        dictionary
    );
}
}
}

Foam::wordList
Foam::radiationModels::absorptionEmissionModels::greyMean::readSpecies(const dictionary& dict)
{
    wordList names;
    forAllConstIter(dictionary, dict, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }
        const word& key = iter().keyword();
        names.append(key);
    }

    return names;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::greyMean::greyMean
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& modelName
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(modelName + "Coeffs")),
    speciesNames_(readSpecies(coeffsDict_)),
    coeffs_(speciesNames_.size(), absorptionCoeffs()),
    eos_(speciesNames_.size()),
    Yj_(speciesNames_.size()),
    p_(mesh.lookupObject<volScalarField>("p")),
    T_(mesh.lookupObject<volScalarField>("T"))
{
    forAll(speciesNames_, i)
    {
        word specieName = speciesNames_[i];
        coeffs_[i].initialise(coeffsDict_.subDict(specieName));
        eos_.set
        (
            i,
            &mesh.lookupObjectRef<fluidThermoModel>
            (
                IOobject::groupName("basicThermo", specieName)
            )
        );
        if (mesh.foundObject<volScalarField>(IOobject::groupName("Y", specieName)))
        {
            Yj_.set
            (
                i,
                &mesh.lookupObjectRef<volScalarField>(IOobject::groupName("Y", specieName))
            );
        }
        else
        {
            Yj_.set
            (
                i,
                &mesh.lookupObjectRef<volScalarField>
                (
                    IOobject::groupName("alpha", specieName)
                )
            );
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::greyMean::~greyMean()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::greyMean::aCont
(
    const label bandI
) const
{
    tmp<volScalarField> ta
    (
        volScalarField::New
        (
            "aCont" + name(bandI),
            mesh(),
            dimensionedScalar(dimless/dimLength, 0),
            extrapolatedCalculatedFvPatchVectorField::typeName
        )
    );

    scalarField& a = ta.ref().primitiveFieldRef();

    tmp<volScalarField> invWt;
    bool vf = Yj_[0].member() == "alpha";
    if (vf)
    {
        invWt = tmp<volScalarField>
        (
            new volScalarField(T_.mesh().lookupObject<volScalarField>("rho"))
        );
    }
    else
    {
        invWt = Yj_[0]/eos_[0].W();
        for (label s = 1; s < Yj_.size(); s++)
        {
            invWt.ref() += Yj_[s]/eos_[s].W();
        }
    }

    forAll(a, celli)
    {
        forAll(Yj_, s)
        {
            scalar Xipi = 0.0;

            scalar Xk;
            if (vf)
            {
                Xk = Yj_[s][celli]*eos_[s].rho()[celli]/invWt()[celli];
            }
            else
            {
                Xk = Yj_[s][celli]/(eos_[s].Wi(celli)*invWt()[celli]);
            }
            Xipi = Xk*paToAtm(p_[celli]);

            const absorptionCoeffs::coeffArray& b = coeffs_[s].coeffs(T_[celli]);

            scalar Ti = max(T_[celli], small);
            // negative temperature exponents
            if (coeffs_[s].invTemp())
            {
                Ti = 1.0/T_[celli];
            }
            a[celli] +=
                Xipi
               *(
                    ((((b[5]*Ti + b[4])*Ti + b[3])*Ti + b[2])*Ti + b[1])*Ti
                  + b[0]
                );
        }
    }
    a = max(a, small);
    ta.ref().correctBoundaryConditions();
    return ta;
}


Foam::scalar
Foam::radiationModels::absorptionEmissionModels::greyMean::aConti
(
    const label celli,
    const label bandI
) const
{
    scalar a = 0.0;

    scalar invWt = Yj_[0][celli]/eos_[0].Wi(celli);
    for (label s = 1; s < Yj_.size(); s++)
    {
        invWt += Yj_[s][celli]/eos_[s].Wi(celli);
    }

    forAll(Yj_, s)
    {
        scalar Xk = Yj_[s][celli]/(eos_[s].Wi(celli)*invWt);
        scalar Xipi = Xk*paToAtm(p_[celli]);

        const absorptionCoeffs::coeffArray& b = coeffs_[s].coeffs(T_[celli]);

        scalar Ti = max(T_[celli], small);
        // negative temperature exponents
        if (coeffs_[s].invTemp())
        {
            Ti = 1.0/T_[celli];
        }
        a +=
            Xipi
            *(
                ((((b[5]*Ti + b[4])*Ti + b[3])*Ti + b[2])*Ti + b[1])*Ti
                + b[0]
            );
    }
    a = max(a, small);
    return a;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::greyMean::eCont
(
    const label bandI
) const
{
    return aCont(bandI);
}


Foam::scalar
Foam::radiationModels::absorptionEmissionModels::greyMean::eConti
(
    const label celli,
    const label bandI
) const
{
    return aConti(celli, bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::greyMean::ECont
(
    const label bandI
) const
{
    return absorptionEmissionModel::ECont(bandI);
}


Foam::scalar
Foam::radiationModels::absorptionEmissionModels::greyMean::EConti
(
    const label celli,
    const label bandI
) const
{
    return absorptionEmissionModel::EConti(celli, bandI);
}


// ************************************************************************* //
