/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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

#define checkTimeIndex_ checkTimeIndex_; public:
#include "fvConstraints.H"
#undef checkTimeIndex_

#include "IBMFvModel.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(IBMForce, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        IBMForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::IBMForce::readCoeffs()
{
    phaseNames_ =
        fvModel::coeffs().lookupOrDefault<wordList>
        (
            "phases",
            wordList(1, word::null)
        );
}


Foam::tmp<Foam::volScalarField>
Foam::fv::IBMForce::alphaRho(const word& phase) const
{
    if (phase == word::null)
    {
        // Assuming incompressible
        if (!mesh().foundObject<volScalarField>("rho"))
        {
            return volScalarField::New
            (
                "rho",
                mesh(),
                dimensionedScalar(dimless, 1.0)
            );
        }

        return mesh().lookupObject<volScalarField>("rho");
    }
    else
    {
        word alphaName(IOobject::groupName("alpha", phase));
        word rhoName(IOobject::groupName("rho", phase));
        return
            mesh().lookupObject<volScalarField>(alphaName)
          * mesh().lookupObject<volScalarField>(rhoName);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::fv::IBMForce::alphaRhoOld(const word& phase) const
{
    if (phase == word::null)
    {
        // Assuming incompressible
        if (mesh().foundObject<volScalarField>("rho"))
        {
            return volScalarField::New
            (
                "rho",
                mesh(),
                dimensionedScalar(dimless, 1.0)
            );
        }

        return mesh().lookupObject<volScalarField>("rho").oldTime();
    }
    else
    {
        word alphaName(IOobject::groupName("alpha", phase));
        word rhoName(IOobject::groupName("rho", phase));
        return
            mesh().lookupObject<volScalarField>(alphaName).oldTime()
          * mesh().lookupObject<volScalarField>(rhoName).oldTime();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::IBMForce::IBMForce
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    fvConstraint(name, modelType, dict, mesh),
    phaseNames_(),
    ibm_(mesh)
{
    readCoeffs();
    fvConstraints& constraints =
        fvConstraints::New(const_cast<fvMesh&>(mesh));
    dictionary& constraintDict(constraints);
    constraintDict.add(name, fvModel::coeffs());
    PtrListDictionary<fvConstraint>& constraintList(constraints);
    constraints.constrainedFields_.setSize(constraintList.size() + 1);
    constraintList.setSize(constraintList.size() + 1);
    constraintList.set
    (
        constraintList.size() - 1,
        name,
        this
    );
    constraints.constrainedFields_.set
    (
        constraintList.size() - 1,
        new wordHashSet()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::IBMForce::constrainedFields() const
{
    wordList fields;
    forAll(phaseNames_, i)
    {
        const word& phase = phaseNames_[i];
        if (ibm_.thermalForcingNeeded())
        {
            const basicThermo& thermo =
                mesh().lookupObject<basicThermo>
                (
                    IOobject::groupName(basicThermo::dictName, phase)
                );
            fields.append(thermo.he().name());
        }
        fields.append(IOobject::groupName("U", phase));
    }
    return fields;
}

// FOR_ALL_FIELD_TYPES
// (
//     IMPLEMENT_FV_CONSTRAINT_CONSTRAIN,
//     fv::IBMForce
// );
IMPLEMENT_FV_CONSTRAINT_CONSTRAIN(scalar, fv::IBMForce);
IMPLEMENT_FV_CONSTRAINT_CONSTRAIN(vector, fv::IBMForce);

void Foam::fv::IBMForce::correct()
{
    ibm_.solve();
}


bool Foam::fv::IBMForce::read(const dictionary& dict)
{
    if (fvConstraint::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
