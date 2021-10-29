/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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

#include "immersedBoundaryObjectListSolver.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::immersedBoundaryObjectListSolver::forcing
(
    const GeometricField<Type, fvPatchField, volMesh>& F,
    const volScalarField& alphaRho,
    const GeometricField<Type, fvPatchField, volMesh>& alphaRhoFOld,
    const GeometricField<Type, fvPatchField, volMesh>& RHS,
    const dimensionedScalar& dt
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tmpForce
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "IBM:Forcing" + F.name(),
                F.time().timeName(),
                F.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            F.mesh(),
            dimensioned<Type>(RHS.dimensions(), Zero)
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& force = tmpForce.ref();
    forAll(objects_, i)
    {
        objects_[i].addForcing
        (
            F.name(), force, alphaRho, alphaRhoFOld, RHS, dt.value()
        );
    }
    return tmpForce;
}


template<class Type>
void Foam::immersedBoundaryObjectListSolver::add
(
    GeometricField<Type, fvPatchField, volMesh>& f
)
{
    const dictionary& defaultDict
    (
        ibmDict_.subDict("defaultBoundaries")
    );

    forAll(objects_, i)
    {
        objects_[i].addField(f, defaultDict);
    }
}


// ************************************************************************* //
