/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::cellCellStencil::createField
(
    const fvMesh& mesh,
    const word& name,
    const UList<Type>& psi
)
{
    tmp<volScalarField> tfld
    (
        new volScalarField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimless, Zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& fld = tfld.ref();

    forAll(psi, cellI)
    {
        fld[cellI] = psi[cellI];
    }
    return tfld;
}


// ************************************************************************* //
