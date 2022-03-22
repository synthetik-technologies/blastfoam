/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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

Class
    Foam::interpolationWeight1D
    Foam::interpolationWeights1D::floor
    Foam::interpolationWeights1D::ceil
    Foam::interpolationWeights1D::linearClamp
    Foam::interpolationWeights1D::linearExtrapolated
    Foam::interpolationWeights1D::quadraticClamp
    Foam::interpolationWeights1D::quadraticExtrapolated
    Foam::interpolationWeights1D::cubicClamp
    Foam::interpolationWeights1D::cubicExtrapolated


Description
    Collection of methods for creating 1D interpolation weights

SourceFiles
    interpolationWeights1D.C

\*---------------------------------------------------------------------------*/

#ifndef interpolationWeights1D_H
#define interpolationWeights1D_H

#include "scalar.H"
#include "List.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class interpolationWeight1D Declaration
\*---------------------------------------------------------------------------*/


class interpolationWeight1D
{
public:

    virtual word type() const = 0;

    virtual ~interpolationWeight1D()
    {}

    static autoPtr<interpolationWeight1D> New(const word&, const label);
    virtual autoPtr<interpolationWeight1D> clone() const = 0;

    static scalar linearWeight
    (
        const scalar x,
        const scalar x0,
        const scalar x1
    )
    {
        return (x - x0)/(x1 - x0);
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        const List<scalar>& xs,
        List<label>& indices,
        List<scalar>& weights
    ) const = 0;
};


namespace interpolationWeights1D
{

/*---------------------------------------------------------------------------*\
                        Floor interpolation Declaration
\*---------------------------------------------------------------------------*/

class floor
:
    public interpolationWeight1D
{
public:

    virtual ~floor()
    {}

    virtual word type() const
    {
        return "floor";
    }

    virtual autoPtr<interpolationWeight1D> clone() const
    {
        return autoPtr<interpolationWeight1D>(new floor());
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        const List<scalar>& xs,
        List<label>& indices,
        List<scalar>& weights
    ) const;
};


/*---------------------------------------------------------------------------*\
                        Ceil interpolation Declaration
\*---------------------------------------------------------------------------*/

class ceil
:
    public interpolationWeight1D
{
public:

    virtual ~ceil()
    {}

    virtual word type() const
    {
        return "ceil";
    }

    virtual autoPtr<interpolationWeight1D> clone() const
    {
        return autoPtr<interpolationWeight1D>(new ceil());
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        const List<scalar>& xs,
        List<label>& indices,
        List<scalar>& weights
    ) const;
};


/*---------------------------------------------------------------------------*\
                        linear interpolation Declaration
\*---------------------------------------------------------------------------*/\

class linearExtrapolated
:
    public interpolationWeight1D
{
public:

    virtual ~linearExtrapolated()
    {}

    virtual word type() const
    {
        return "linearExtrapolated";
    }

    virtual autoPtr<interpolationWeight1D> clone() const
    {
        return autoPtr<interpolationWeight1D>(new linearExtrapolated());
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        const List<scalar>& xs,
        List<label>& indices,
        List<scalar>& weights
    ) const;
};


class linearClamp
:
    public linearExtrapolated
{
public:

    virtual ~linearClamp()
    {}

    virtual word type() const
    {
        return "linearClamp";
    }

    virtual autoPtr<interpolationWeight1D> clone() const
    {
        return autoPtr<interpolationWeight1D>(new linearClamp());
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        const List<scalar>& xs,
        List<label>& indices,
        List<scalar>& weights
    ) const;
};

/*---------------------------------------------------------------------------*\
                        Quadratic interpolation Declaration
\*---------------------------------------------------------------------------*/


class quadraticExtrapolated
:
    public interpolationWeight1D
{
public:

    virtual ~quadraticExtrapolated()
    {}

    virtual word type() const
    {
        return "quadraticExtrapolated";
    }

    virtual autoPtr<interpolationWeight1D> clone() const
    {
        return autoPtr<interpolationWeight1D>(new quadraticExtrapolated());
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        const List<scalar>& xs,
        List<label>& indices,
        List<scalar>& weights
    ) const;
};


class quadraticClamp
:
    public quadraticExtrapolated
{
public:

    virtual ~quadraticClamp()
    {}

    virtual word type() const
    {
        return "quadraticClamp";
    }

    virtual autoPtr<interpolationWeight1D> clone() const
    {
        return autoPtr<interpolationWeight1D>(new quadraticClamp());
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        const List<scalar>& xs,
        List<label>& indices,
        List<scalar>& weights
    ) const;
};


/*---------------------------------------------------------------------------*\
                        Cubic interpolation Declaration
\*---------------------------------------------------------------------------*/


class cubicExtrapolated
:
    public interpolationWeight1D
{
public:

    virtual ~cubicExtrapolated()
    {}

    virtual word type() const
    {
        return "cubicExtrapolated";
    }

    virtual autoPtr<interpolationWeight1D> clone() const
    {
        return autoPtr<interpolationWeight1D>(new cubicExtrapolated());
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        const List<scalar>& xs,
        List<label>& indices,
        List<scalar>& weights
    ) const;
};


class cubicClamp
:
    public cubicExtrapolated
{
public:

    virtual ~cubicClamp()
    {}

    virtual word type() const
    {
        return "cubicClamp";
    }

    virtual autoPtr<interpolationWeight1D> clone() const
    {
        return autoPtr<interpolationWeight1D>(new cubicClamp());
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        const List<scalar>& xs,
        List<label>& indices,
        List<scalar>& weights
    ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interpolationWeights1D


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
