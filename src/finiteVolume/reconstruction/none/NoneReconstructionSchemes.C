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

#include "NoneReconstructionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineNamedTemplateTypeNameAndDebug(NoneReconstructionScheme<scalar>, 0);
ReconstructionScheme<scalar>::adddictionaryConstructorToTable
    <
        NoneReconstructionScheme<scalar>
    > addnoneMUSCLscalardictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug(NoneReconstructionScheme<vector>, 0);
ReconstructionScheme<vector>::adddictionaryConstructorToTable
    <
        NoneReconstructionScheme<vector>
    > addnoneMUSCLvectordictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug(NoneReconstructionScheme<symmTensor>, 0);
ReconstructionScheme<symmTensor>::adddictionaryConstructorToTable
    <
        NoneReconstructionScheme<symmTensor>
    > addnoneMUSCLsymmTensordictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug
(
    NoneReconstructionScheme<sphericalTensor>,
    0
);
ReconstructionScheme<sphericalTensor>::adddictionaryConstructorToTable
    <
        NoneReconstructionScheme<sphericalTensor>
    > addnoneMUSCLsphericalTensordictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug(NoneReconstructionScheme<tensor>, 0);
ReconstructionScheme<tensor>::adddictionaryConstructorToTable
    <
        NoneReconstructionScheme<tensor>
    > addnoneMUSCLtensordictinoaryConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
