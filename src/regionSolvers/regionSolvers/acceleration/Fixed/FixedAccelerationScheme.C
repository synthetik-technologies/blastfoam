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

#include "FixedAccelerationScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::accelerationSchemes::Fixed<Type, Patch, Mesh>::Fixed
(
    GeometricField<Type, Patch, Mesh>& field,
    autoPtr<PatchFieldSelector<Type>> selector,
    const dictionary& dict
)
:
    AccelerationSchemeBase<Type, Patch, Mesh>(typeName, field, selector, dict),
    relaxationFactor_(1.0)
{
    read(dict);
    Info<< indent << "relaxationFactor: " << relaxationFactor_ << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::accelerationSchemes::Fixed<Type, Patch, Mesh>::~Fixed()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
void Foam::accelerationSchemes::Fixed<Type, Patch, Mesh>::relax
(
    const label iter
)
{
    this->updateError();

    Patch<Type>& pf = this->field_.boundaryFieldRef()[this->patchi_];
    Field<Type>& pfield =this->selector_->relaxField(pf);
    const Field<Type>& pfieldPrev =
        this->selector_->relaxField
        (
            this->field_.prevIter().boundaryField()[this->patchi_]
        );
    pfield =
        pfieldPrev*(1.0 - relaxationFactor_) + relaxationFactor_*pfield;
    this->selector_->correct(pf);
}


template<class Type, template<class> class Patch, class Mesh>
void Foam::accelerationSchemes::Fixed<Type, Patch, Mesh>::read
(
    const dictionary& dict
)
{
    this->coeffDict(dict).lookup("relaxationFactor") >> relaxationFactor_;
}


// ************************************************************************* //
