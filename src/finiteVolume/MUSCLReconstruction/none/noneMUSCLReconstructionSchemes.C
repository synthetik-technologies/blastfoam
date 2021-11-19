/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020
     \\/     M anipulation  | Synthetik Applied Technology
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "noneMUSCLReconstructionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineNamedTemplateTypeNameAndDebug(noneMUSCLReconstructionScheme<scalar>, 0);
MUSCLReconstructionScheme<scalar>::adddictionaryConstructorToTable
    <
        noneMUSCLReconstructionScheme<scalar>
    > addnoneMUSCLscalardictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug(noneMUSCLReconstructionScheme<vector>, 0);
MUSCLReconstructionScheme<vector>::adddictionaryConstructorToTable
    <
        noneMUSCLReconstructionScheme<vector>
    > addnoneMUSCLvectordictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug(noneMUSCLReconstructionScheme<symmTensor>, 0);
MUSCLReconstructionScheme<symmTensor>::adddictionaryConstructorToTable
    <
        noneMUSCLReconstructionScheme<symmTensor>
    > addnoneMUSCLsymmTensordictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug
(
    noneMUSCLReconstructionScheme<sphericalTensor>,
    0
);
MUSCLReconstructionScheme<sphericalTensor>::adddictionaryConstructorToTable
    <
        noneMUSCLReconstructionScheme<sphericalTensor>
    > addnoneMUSCLsphericalTensordictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug(noneMUSCLReconstructionScheme<tensor>, 0);
MUSCLReconstructionScheme<tensor>::adddictionaryConstructorToTable
    <
        noneMUSCLReconstructionScheme<tensor>
    > addnoneMUSCLtensordictinoaryConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
