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

#include "upwindMUSCLReconstructionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineNamedTemplateTypeNameAndDebug(upwindMUSCLReconstructionScheme<scalar>, 0);
MUSCLReconstructionScheme<scalar>::adddictionaryConstructorToTable
    <
        upwindMUSCLReconstructionScheme<scalar>
    > addupwindMUSCLscalardictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug(upwindMUSCLReconstructionScheme<vector>, 0);
MUSCLReconstructionScheme<vector>::adddictionaryConstructorToTable
    <
        upwindMUSCLReconstructionScheme<vector>
    > addupwindMUSCLvectordictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug(upwindMUSCLReconstructionScheme<symmTensor>, 0);
MUSCLReconstructionScheme<symmTensor>::adddictionaryConstructorToTable
    <
        upwindMUSCLReconstructionScheme<symmTensor>
    > addupwindMUSCLsymmTensordictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug
(
    upwindMUSCLReconstructionScheme<sphericalTensor>,
    0
);
MUSCLReconstructionScheme<sphericalTensor>::adddictionaryConstructorToTable
    <
        upwindMUSCLReconstructionScheme<sphericalTensor>
    > addupwindMUSCLsphericalTensordictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug(upwindMUSCLReconstructionScheme<tensor>, 0);
MUSCLReconstructionScheme<tensor>::adddictionaryConstructorToTable
    <
        upwindMUSCLReconstructionScheme<tensor>
    > addupwindMUSCLtensordictinoaryConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
