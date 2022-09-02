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

#include "UpwindMUSCLReconstructionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineNamedTemplateTypeNameAndDebug(UpwindMUSCLReconstructionScheme<scalar>, 0);
ReconstructionScheme<scalar>::adddictionaryConstructorToTable
    <
        UpwindMUSCLReconstructionScheme<scalar>
    > addupwindMUSCLscalardictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug(UpwindMUSCLReconstructionScheme<vector>, 0);
ReconstructionScheme<vector>::adddictionaryConstructorToTable
    <
        UpwindMUSCLReconstructionScheme<vector>
    > addupwindMUSCLvectordictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug(UpwindMUSCLReconstructionScheme<symmTensor>, 0);
ReconstructionScheme<symmTensor>::adddictionaryConstructorToTable
    <
        UpwindMUSCLReconstructionScheme<symmTensor>
    > addupwindMUSCLsymmTensordictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug
(
    UpwindMUSCLReconstructionScheme<sphericalTensor>,
    0
);
ReconstructionScheme<sphericalTensor>::adddictionaryConstructorToTable
    <
        UpwindMUSCLReconstructionScheme<sphericalTensor>
    > addupwindMUSCLsphericalTensordictinoaryConstructorToTable_;

defineNamedTemplateTypeNameAndDebug(UpwindMUSCLReconstructionScheme<tensor>, 0);
ReconstructionScheme<tensor>::adddictionaryConstructorToTable
    <
        UpwindMUSCLReconstructionScheme<tensor>
    > addupwindMUSCLtensordictinoaryConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
