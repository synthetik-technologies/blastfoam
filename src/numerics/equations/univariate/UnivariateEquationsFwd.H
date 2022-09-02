/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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

#ifndef UnivariateEquationsFwd_H
#define UnivariateEquationsFwd_H

#include "UnivariateEquation.H"
#include "adaptiveTypesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef univariateEquation<scalar> scalarUnivariateEquation;
typedef univariateEquation<vector> vectorUnivariateEquation;
typedef univariateEquation<symmTensor> symmTensorUnivariateEquation;
typedef univariateEquation<sphericalTensor>
    sphericalTensorUnivariateEquation;
typedef univariateEquation<tensor> tensorUnivariateEquation;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef UnivariateEquation<scalar> ScalarUnivariateEquation;
typedef UnivariateEquation<vector> VectorUnivariateEquation;
typedef UnivariateEquation<symmTensor> SymmTensorUnivariateEquation;
typedef UnivariateEquation<sphericalTensor>
    SphericalTensorUnivariateEquation;
typedef UnivariateEquation<tensor> TensorUnivariateEquation;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
