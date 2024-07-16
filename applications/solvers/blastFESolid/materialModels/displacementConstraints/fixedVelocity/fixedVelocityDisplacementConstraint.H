/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
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

Class
    Foam::displacementConstraints::fixedVelocity

Description
    Fixed velocity nodal constraint

SourceFiles
    fixedVelocityNodeConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef fixedVelocityNodeConstraint_H
#define fixedVelocityNodeConstraint_H

#include "displacementConstraint.H"
#include "Function1.H"
#include "pointFields.H"

namespace Foam
{
namespace displacementConstraints
{

class fixedVelocity
:
    public displacementConstraint
{
protected:

    //- Displacement function
    autoPtr<Function1<vector>> velocity_;

    fixedVelocity(const fixedVelocity&);

public:

    TypeName("fixedVelocity");

    fixedVelocity
    (
        const word& name,
        pointVectorField& D,
        pointVectorField& U,
        const dictionary& dict
    );

    virtual autoPtr<displacementConstraint> clone() const
    {
        return autoPtr<displacementConstraint>
        (
            new fixedVelocity(*this)
        );
    }

    virtual ~fixedVelocity();

    // Evaluate the coefficient at the given point
    virtual void constrain();
};

} // End namespace displacementConstraints
} // End namespace Foam

#endif
