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

#include "QNBaseAccelerationScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::accelerationSchemes::QNBase<Type, Patch, Mesh>::QNBase
(
    const word& type,
    GeometricField<Type, Patch, Mesh>& field,
    autoPtr<PatchFieldSelector<Type>> selector,
    const dictionary& dict
)
:
    ResidualBase<Type, Patch, Mesh>(type, field, selector, dict),

    nCouplingTimes_(0),

    Vs_(),
    Ws_(),
    times_()
{
    read(dict);
    Info<< indent << "nCouplingTimes: " << nCouplingTimes_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::accelerationSchemes::QNBase<Type, Patch, Mesh>::~QNBase()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
void Foam::accelerationSchemes::QNBase<Type, Patch, Mesh>::updateVW
(
    const label iter
)
{
    const Field<Type>& pfield =
        this->selector_->relaxField
        (
            this->field_.boundaryField()[this->patchi_]
        );

    if (iter == 0)
    {
    }
    else if (iter > 0)
    {
        Vs_.append(this->residuals_ - this->VRef());
        Ws_.append(pfield - this->WRef());
        times_.append(this->field_.time().timeIndex());
    }
}


template<class Type, template<class> class Patch, class Mesh>
void Foam::accelerationSchemes::QNBase<Type, Patch, Mesh>::clear
(
    const bool full
)
{
    ResidualBase<Type, Patch, Mesh>::clear(full);

    if (!nCouplingTimes_ | full)
    {
        Vs_.clear();
        Ws_.clear();
        times_.clear();

        return;
    }

    label startI = times_.size();
    forAll(times_, ti)
    {
        if ((this->field_.time().timeIndex() - times_[ti]) < nCouplingTimes_)
        {
            startI = ti;
            break;
        }
    }

    const label oldSize = times_.size();
    if (startI != times_.size() && startI)
    {
        for (label ti = 0; ti < oldSize-startI; ti++)
        {
            Vs_[ti].transfer(Vs_[ti+startI]);
            Ws_[ti].transfer(Ws_[ti+startI]);
            times_[ti] = times_[ti+startI];
        }
        Vs_.setSize(oldSize-startI);
        Ws_.setSize(oldSize-startI);
        times_.setSize(oldSize-startI);
    }
}


template<class Type, template<class> class Patch, class Mesh>
void Foam::accelerationSchemes::QNBase<Type, Patch, Mesh>::read
(
    const dictionary& dict
)
{
    ResidualBase<Type, Patch, Mesh>::read(dict);
    this->coeffDict(dict).lookup("nCouplingTimes") >> nCouplingTimes_;
}


// ************************************************************************* //
