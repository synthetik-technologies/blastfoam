/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "LinearGeomSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IncrementalModel>
LinearGeomSolid<IncrementalModel>::LinearGeomSolid
(
    const word& type,
    dynamicFvMesh& mesh,
    const bool isSolid
)
:
    IncrementalModel(type, mesh, nonLinGeom(), isSolid)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IncrementalModel>
void LinearGeomSolid<IncrementalModel>::update(const bool correctSigma)
{
    IncrementalModel::updateDisplacement();

    if (correctSigma)
    {
        this->mechanical().correct(this->sigma());
        this->sigma().correctBoundaryConditions();
    }
}


template<class IncrementalModel>
bool LinearGeomSolid<IncrementalModel>::write(const bool write) const
{
    bool good = IncrementalModel::write(write);
    if (write)
    {
        // Write strain fields

        // Total strain
        volSymmTensorField epsilon("epsilon", symm(this->gradD()));

        // Equivalent strain
        volScalarField epsilonEq
        (
            "epsilonEq", sqrt((2.0/3.0)*magSqr(dev(epsilon)))
        );
        Info<< "Max epsilonEq = " << gMax(epsilonEq) << endl;

        good = epsilon.write() && epsilonEq.write();
    }
    return good && IncrementalModel::write(write);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels
} // End namespace Foam

// ************************************************************************* //
