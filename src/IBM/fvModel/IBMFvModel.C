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
#include "IBMFvConstraint.H"

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
    // Check if the fvConstraint exists for the IBM forcing and
    // add it if it does not
    fvConstraints& constraints =
        fvConstraints::New(const_cast<fvMesh&>(mesh));
    dictionary& constraintDict(constraints);
    if (!constraintDict.found(name))
    {
        constraintDict.add(name, coeffs());
        PtrListDictionary<fvConstraint>& constraintList(constraints);
        const label ci = constraintList.size();
        constraints.constrainedFields_.setSize(ci + 1);
        constraintList.setSize(ci + 1);
        constraintList.set
        (
            ci,
            name,
            new IBMForceConstraint
            (
                name,
                modelType,
                coeffs(),
                mesh
            )
        );
        constraints.constrainedFields_.set
        (
            ci,
            new wordHashSet()
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::fv::IBMForceModel::correct()
{
    ibm_.solve();
}


// ************************************************************************* //
