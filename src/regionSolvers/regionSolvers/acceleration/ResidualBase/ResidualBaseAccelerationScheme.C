/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "ResidualBaseAccelerationScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::accelerationSchemes::ResidualBase<Type, Patch, Mesh>::ResidualBase
(
    const word& type,
    GeometricField<Type, Patch, Mesh>& field,
    autoPtr<PatchFieldSelector<Type>> selector,
    const dictionary& dict
)
:
    AccelerationSchemeBase<Type, Patch, Mesh>(type, field, selector, dict),

    residuals_(),
    prevResiduals_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::accelerationSchemes::ResidualBase<Type, Patch, Mesh>::~ResidualBase()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
void Foam::accelerationSchemes::ResidualBase<Type, Patch, Mesh>::updateError()
{
    this->error_ = Zero;

    const Field<Type>& pfield = this->selector_->field();
    const Field<Type>& pfieldPrev =
        this->selector_->relaxField
        (
            this->field_.prevIter().boundaryField()[this->patchi_]
        );

    prevResiduals_.transfer(residuals_);
    residuals_ = pfield - pfieldPrev;

    scalar maxErrorMagSqr = 0.0;
    const label n = returnReduce(pfield.size(), sumOp<label>());
    forAll(residuals_, i)
    {
        scalar errorMagSqr = magSqr(residuals_[i]);
        this->error_ += magSqr(errorMagSqr);
        maxErrorMagSqr = max(maxErrorMagSqr, errorMagSqr);
    }

    reduce(this->error_, sumOp<scalar>());
    reduce(maxErrorMagSqr, maxOp<scalar>());
    this->error_ = sqrt(this->error_/(maxErrorMagSqr + small))/scalar(n);

    if (this->initialError_ < 0)
    {
        this->initialError_ = this->error_;
    }
}


template<class Type, template<class> class Patch, class Mesh>
void Foam::accelerationSchemes::ResidualBase<Type, Patch, Mesh>::clear
(
    const bool full
)
{
    AccelerationSchemeBase<Type, Patch, Mesh>::clear(full);
    residuals_ = Zero;
}

// ************************************************************************* //
