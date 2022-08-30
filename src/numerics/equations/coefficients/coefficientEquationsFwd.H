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

#ifndef coefficientEquationsFwd_H
#define coefficientEquationsFwd_H

#include "Coefficients.H"
#include "coefficientsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef CoefficientEquation<scalar> ScalarCoefficientEquation;
typedef CoefficientEquation<vector> VectorCoefficientEquation;
typedef CoefficientEquation<symmTensor> SymmTensorCoefficientEquation;
typedef CoefficientEquation<sphericalTensor> SphericalTensorCoefficientEquation;
typedef CoefficientEquation<tensor> TensorCoefficientEquation;


typedef LinearEquation<scalar> ScalarLinearEquation;
typedef LinearEquation<vector> VectorLinearEquation;
typedef LinearEquation<symmTensor> SymmTensorLinearEquation;
typedef LinearEquation<sphericalTensor> SphericalTensorLinearEquation;
typedef LinearEquation<tensor> TensorLinearEquation;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef UnivariateCoefficientEquation<scalar>
    ScalarUnivariateCoefficientEquation;
typedef UnivariateCoefficientEquation<vector>
    VectorUnivariateCoefficientEquation;
typedef UnivariateCoefficientEquation<symmTensor>
    SymmTensorUnivariateCoefficientEquation;
typedef UnivariateCoefficientEquation<sphericalTensor>
    SphericalTensorUnivariateCoefficientEquation;
typedef UnivariateCoefficientEquation<tensor>
    TensorUnivariateCoefficientEquation;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef LinearUnivariateEquation<scalar> ScalarLinearUnivariateEquation;
typedef LinearUnivariateEquation<vector> VectorLinearUnivariateEquation;
typedef LinearUnivariateEquation<symmTensor> SymmTensorLinearUnivariateEquation;
typedef LinearUnivariateEquation<sphericalTensor>
    SphericalTensorLinearUnivariateEquation;
typedef LinearUnivariateEquation<tensor> TensorLinearUnivariatEquation;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
