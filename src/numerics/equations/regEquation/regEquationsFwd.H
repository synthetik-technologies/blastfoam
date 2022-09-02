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

#ifndef regEquationsFwd_H
#define regEquationsFwd_H

#include "regEquation.H"
#include "EquationsFwd.H"
#include "UnivariateEquationsFwd.H"
#include "MultivariateEquationsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef regEquation<scalar, Equation> ScalarRegEquation;
typedef regEquation<vector, Equation> VectorRegEquation;
typedef regEquation<symmTensor, Equation> SymmTensorRegEquation;
typedef regEquation<sphericalTensor, Equation> SphericalTensorRegEquation;
typedef regEquation<tensor, Equation> TensorRegEquation;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef regEquation<scalar, UnivariateEquation>
    ScalarRegUnivariateEquation;
typedef regEquation<vector, UnivariateEquation>
    VectorRegUnivariateEquation;
typedef regEquation<symmTensor, UnivariateEquation>
    SymmTensorRegUnivariateEquation;
typedef regEquation<sphericalTensor, UnivariateEquation>
    SphericalTensorRegUnivariateEquation;
typedef regEquation<tensor, UnivariateEquation>
    TensorRegUnivariateEquation;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef regEquation<scalar, MultivariateEquation>
    ScalarRegMultivariateEquation;
typedef regEquation<vector, MultivariateEquation>
    VectorRegMultivariateEquation;
typedef regEquation<symmTensor, MultivariateEquation>
    SymmTensorRegMultivariateEquation;
typedef regEquation<sphericalTensor, MultivariateEquation>
    SphericalTensorRegMultivariateEquation;
typedef regEquation<tensor, MultivariateEquation>
    TensorRegMultivariateEquation;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
