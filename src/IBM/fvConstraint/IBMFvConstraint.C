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
#include "fvModels.H"
#undef checkTimeIndex_

#include "IBMFvConstraint.H"
#include "IBMFvModel.H"

#include "fvMatrices.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(IBMForceConstraint, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::IBMForceConstraint::readCoeffs()
{
    forcedFields_ = coeffs().lookupOrDefault<wordList>
    (
        "forcedFields",
        wordList()
    );
}


Foam::tmp<Foam::volScalarField>
Foam::fv::IBMForceConstraint::alphaRho(const word& phase) const
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::IBMForceConstraint::IBMForceConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvConstraint(name, modelType, dict, mesh),
    ibm_
    (
        immersedBoundaryObjectListSolver::New
        (
            const_cast<fvMesh&>(mesh)
        )
    ),
    forcedFields_()
{
    readCoeffs();
    fvModels& models =
        fvModels::New(const_cast<fvMesh&>(mesh));
    dictionary& modelDict(models);
    if (!modelDict.found(name))
    {
        modelDict.add(name, coeffs());
        PtrListDictionary<fvModel>& modelList(models);
        label mi = modelList.size();
        models.addSupFields_.setSize(mi + 1);
        modelList.setSize(mi + 1);
        modelList.set
        (
            mi,
            name,
            new IBMForceModel
            (
                name,
                modelType,
                coeffs(),
                mesh
            )
        );
        models.addSupFields_.set(mi, new wordHashSet());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList
Foam::fv::IBMForceConstraint::constrainedFields() const
{
    return forcedFields_;
}

FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_CONSTRAINT_CONSTRAIN,
    fv::IBMForceConstraint
);


bool Foam::fv::IBMForceConstraint::read(const dictionary& dict)
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
