/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "blastGreyMean.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "multicomponentThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace absorptionEmissionModels
{
    defineTypeNameAndDebug(blastGreyMean, 0);

    addToRunTimeSelectionTable
    (
        blastAbsorptionEmissionModel,
        blastGreyMean,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::blastGreyMean::blastGreyMean
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& modelName
)
:
    blastAbsorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(modelName + "Coeffs")),
    speciesNames_(0),
    specieIndex_(),
    lookUpTablePtr_(),
    thermo_(mesh.lookupObject<fluidThermo>(physicalProperties::typeName)),
    Yj_()
{
    label nFunc = 0;
    label nPhases = 0;
    HashTable<label> phases;
    bool anyPhase = false;
    bool noPhase = false;
    forAllConstIter(dictionary, coeffsDict_, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }
        const word& key = iter().keyword();
        const dictionary& dict = iter().dict();
        if (!speciesNames_.insert(key, nFunc))
        {
            FatalIOErrorInFunction(dict)
                << "Multiple entries named " << key << "exist"
                << abort(FatalIOError);
        }
        coeffs_.append(absorptionCoeffs());
        coeffs_.last().initialise(dict);
        if (dict.found("phase"))
        {
            const word phaseName(dict.lookup("phase"));
            if (!phases.found(phaseName))
            {
                phases.insert(phaseName, nPhases++);
            }
            speciePhases_.append({phases[phaseName]});
            anyPhase = true;
        }
        else if (dict.found("phases"))
        {
            const wordList phaseNames(dict.lookup("phases"));
            speciePhases_.setSize(speciePhases_.size() + 1);
            forAll(phaseNames, i)
            {
                const word& phaseName = phaseNames[i];
                if (!phases.found(phaseName))
                {
                    phases.insert(phaseName, nPhases++);
                }
                speciePhases_.last().append(phases[phaseName]);
            }
            anyPhase = true;
        }
        else
        {
            noPhase = true;
        }
        nFunc++;
    }
    if (anyPhase && noPhase)
    {
        FatalErrorInFunction
            << "Some species have a phase specified, and some do not."
            << " Either all species should have a phase or none." << endl
            << abort(FatalError);
    }
    else if (anyPhase)
    {
        volumeFractions_.setSize(phases.size());
        forAllConstIter(HashTable<label>, phases, iter)
        {
            volumeFractions_.set
            (
                iter(),
                &mesh.lookupObjectRef<volScalarField>
                (
                    IOobject::groupName("alpha", iter.key())
                )
            );
        }
    }

    if (isA<multicomponentThermo>(thermo_))
    {
        mixture_.set(dynamic_cast<const multicomponentThermo*>(&thermo_));

        if (coeffsDict_.found("lookUpTableFileName"))
        {
            const word name = coeffsDict_.lookup("lookUpTableFileName");
            if (name != "none")
            {
                lookUpTablePtr_.set
                (
                    new interpolationLookUpTable
                    (
                        fileName(coeffsDict_.lookup("lookUpTableFileName")),
                        mesh.time().constant(),
                        mesh
                    )
                );

                if (!mesh.foundObject<volScalarField>("ft"))
                {
                    FatalErrorInFunction
                        << "specie ft is not present to use with "
                        << "lookUpTableFileName " << nl
                        << exit(FatalError);
                }
                ft_.set(&mesh.lookupObject<volScalarField>("ft"));
            }
        }

        // Check that all the species on the dictionary are present in the
        // look-up table and save the corresponding indices of the look-up table

        label j = 0;
        forAllConstIter(HashTable<label>, speciesNames_, iter)
        {
            if (!lookUpTablePtr_.empty())
            {
                if (lookUpTablePtr_().found(iter.key()))
                {
                    label index = lookUpTablePtr_().findFieldIndex(iter.key());

                    Info<< "specie: " << iter.key() << " found on look-up table "
                        << " with index: " << index << endl;

                    specieIndex_[iter()] = index;
                }
                else if (mesh.foundObject<volScalarField>(iter.key()))
                {
                    Yj_.set
                    (
                        j,
                        &mesh.lookupObjectRef<volScalarField>(iter.key())
                    );
                    specieIndex_[iter()] = 0;
                    j++;
                    Info<< "specie: " << iter.key()
                        << " is being solved" << endl;
                }
                else
                {
                    FatalErrorInFunction
                        << "specie: " << iter.key()
                        << " is neither in look-up table: "
                        << lookUpTablePtr_().tableName()
                        << " nor is being solved" << nl
                        << exit(FatalError);
                }
            }
            else if (mesh.foundObject<volScalarField>(iter.key()))
            {
                Yj_.set(j, &mesh.lookupObjectRef<volScalarField>(iter.key()));
                specieIndex_[iter()] = 0;
                j++;
            }
            else
            {
                FatalErrorInFunction
                    << " there is not lookup table and the specie" << nl
                    << iter.key() << nl
                    << " is not found " << nl
                    << exit(FatalError);

            }
        }
    }
    else
    {
        Xj_.setSize(speciesNames_.size());
        forAllConstIter(HashTable<label>, speciesNames_, iter)
        {
            const dictionary& dict = coeffsDict_.subDict(iter.key());
            Xj_[iter()] = dict.lookup<scalar>("molFraction");
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::absorptionEmissionModels::blastGreyMean::~blastGreyMean()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::blastGreyMean::aCont
(
    const label bandI
) const
{
    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();


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

    const unitConversion& unitAtm = units()["atm"];

    scalarField& a = ta.ref().primitiveFieldRef();
    forAll(a, celli)
    {
        forAllConstIter(HashTable<label>, speciesNames_, iter)
        {
            label n = iter();
            scalar Xi = 0.0;
            if (Xj_.size())
            {
                Xi = Xj_[n];
            }
            else if (specieIndex_[n] != 0)
            {
                // Specie found in the lookUpTable.
                const List<scalar>& Ynft =
                    lookUpTablePtr_().lookUp(ft_()[celli]);

                // moles x pressure [atm]
                Xi = Ynft[specieIndex_[n]];
            }
            else
            {
                scalar invWt = 0.0;
                forAll(mixture_->Y(), s)
                {
                    invWt += mixture_->Y(s)[celli]/mixture_->WiValue(s);
                }

                label index = mixture_->species()[iter.key()];
                Xi = mixture_->Y(index)[celli]/(mixture_->WiValue(index)*invWt);
            }

            if (speciePhases_[n].size())
            {
                scalar vf = 0.0;
                forAll(speciePhases_[n], phasei)
                {
                    vf += volumeFractions_[speciePhases_[n][phasei]][celli];
                }
                Xi *= vf;
            }


            scalar Ti = T[celli];

            const absorptionCoeffs::coeffArray& b = coeffs_[n].coeffs(Ti);

            // negative temperature exponents
            if (coeffs_[n].invTemp())
            {
                Ti = 1.0/Ti;
            }
            a[celli] +=
                Xi*unitAtm.toUser(p[celli])
               *(
                    ((((b[5]*Ti + b[4])*Ti + b[3])*Ti + b[2])*Ti + b[1])*Ti
                  + b[0]
                );
        }
    }

    ta.ref().correctBoundaryConditions();
    return ta;
}


Foam::scalar
Foam::radiationModels::absorptionEmissionModels::blastGreyMean::cellaCont
(
    const label celli,
    const label bandI
) const
{
    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();

    scalar a = 0.0;

    const unitConversion& unitAtm = units()["atm"];

    forAllConstIter(HashTable<label>, speciesNames_, iter)
    {
        label n = iter();
        scalar Xi = 0.0;
        if (Xj_.size())
        {
            Xi = Xj_[n];
        }
        else if (specieIndex_[n] != 0)
        {
            // Specie found in the lookUpTable.
            const List<scalar>& Ynft = lookUpTablePtr_().lookUp(ft_()[celli]);

            // moles x pressure [atm]
            Xi = Ynft[specieIndex_[n]];
        }
        else
        {
            scalar invWt = 0.0;
            forAll(mixture_->Y(), s)
            {
                invWt += mixture_->Y(s)[celli]/mixture_->WiValue(s);
            }

            label index = mixture_->species()[iter.key()];
            Xi = mixture_->Y(index)[celli]/(mixture_->WiValue(index)*invWt);
        }

        scalar Ti = T[celli];

        const absorptionCoeffs::coeffArray& b = coeffs_[n].coeffs(Ti);

        // negative temperature exponents
        if (coeffs_[n].invTemp())
        {
            Ti = 1.0/Ti;
        }
        a +=
            Xi*unitAtm.toUser(p[celli])
           *(
                ((((b[5]*Ti + b[4])*Ti + b[3])*Ti + b[2])*Ti + b[1])*Ti
              + b[0]
            );
    }
    return a;
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::blastGreyMean::eCont
(
    const label bandI
) const
{
    return aCont(bandI);
}


Foam::scalar
Foam::radiationModels::absorptionEmissionModels::blastGreyMean::celleCont
(
    const label celli,
    const label bandI
) const
{
    return cellaCont(celli, bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::absorptionEmissionModels::blastGreyMean::ECont
(
    const label bandI
) const
{
    return absorptionEmissionModel::ECont(bandI);
}


Foam::scalar
Foam::radiationModels::absorptionEmissionModels::blastGreyMean::cellECont
(
    const label celli,
    const label bandI
) const
{
    return blastAbsorptionEmissionModel::cellECont(celli, bandI);
}


// ************************************************************************* //
