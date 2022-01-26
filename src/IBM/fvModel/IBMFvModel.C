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

#include "IBMFvModel.H"
#include "IBMFvConstraint.H"

#include "fvConstraints.H"

#include "fvMatrices.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(IBMForceModel, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        IBMForceModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::IBMForceModel::IBMForceModel
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(name, modelType, dict, mesh),
    ibm_
    (
        immersedBoundaryObjectListSolver::New
        (
            const_cast<fvMesh&>(mesh)
        )
    )
{
    // Initialize the shapes
    forAll(ibm_.objects(), i)
    {
        ibm_.objects()[i].shape();
    }

    // Create a temporary fvConstraints object to add the IBMConstraint to
    fvConstraints constraints(mesh);
    IOobject constraintsIO(constraints);
    constraintsIO.registerObject() = false;
    if (constraints.headerOk())
    {
        constraintsIO.readOpt() = IOobject::READ_IF_PRESENT;
    }
    else
    {
        constraintsIO.readOpt() = IOobject::NO_READ;
        constraintsIO.instance() = mesh.time().system();
    }
    IOdictionary constraintsDict
    (
        constraintsIO,
        constraints
    );
    {
        constraintsDict.add(name, coeffs());
        constraintsDict.regIOobject::write();
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::fv::IBMForceModel::correct()
{
    ibm_.solve();
}


// ************************************************************************* //
