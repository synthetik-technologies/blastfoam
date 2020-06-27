/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "fvMotionSolver.H"
#include "fixedValuePointPatchFields.H"
#include "cellMotionFvPatchFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::fvMotionSolver::cellMotionBoundaryTypes
(
    const typename GeometricField<Type, pointPatchField, pointMesh>::
    Boundary& pmUbf
) const
{
    wordList cmUbf = pmUbf.types();

    // Remove global patches from the end of the list
    cmUbf.setSize(fvMesh_.boundary().size());

    forAll(cmUbf, patchi)
    {
        if (isA<fixedValuePointPatchField<Type>>(pmUbf[patchi]))
        {
            cmUbf[patchi] = cellMotionFvPatchField<Type>::typeName;
        }

        if (debug)
        {
            Pout<< "Patch:" << fvMesh_.boundary()[patchi].patch().name()
                << " pointType:" << pmUbf.types()[patchi]
                << " cellType:" << cmUbf[patchi] << endl;
        }
    }

    return cmUbf;
}


// ************************************************************************* //
