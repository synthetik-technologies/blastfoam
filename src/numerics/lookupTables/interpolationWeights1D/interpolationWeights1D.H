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
protected:
    const List<scalar>& xs_;

public:

    virtual word type() const = 0;

    interpolationWeight1D(const List<scalar>& xs)
    :
        xs_(xs)
    {}

    virtual ~interpolationWeight1D()
    {}

    static autoPtr<interpolationWeight1D> New
    (
        const word&,
        const List<scalar>&
    );
    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const = 0;

    virtual bool validate(const bool fail = true) const;

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
        List<label>& indices,
        List<scalar>& weights
    ) const = 0;
};


namespace interpolationWeights1D
{

/*---------------------------------------------------------------------------*\
                        Clamp interpolation Declaration
\*---------------------------------------------------------------------------*/

template<class Interp>
class Clamp
:
    public Interp
{
public:

    Clamp(const List<scalar>& xs)
    :
        Interp(xs)
    {}

    virtual ~Clamp()
    {}

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        List<label>& indices,
        List<scalar>& weights
    ) const
    {
        if (x < this->xs_[0])
        {
            indices.setSize(1);
            weights.setSize(1);

            indices[0] = 0;
            weights[0] = 1.0;

            return;
        }
        else if (x > this->xs_.last())
        {
            indices.setSize(1);
            weights.setSize(1);

            indices[0] = this->xs_.size() - 1;
            weights[0] = 1.0;
            return;
        }
        Interp::updateWeights(x, i, indices, weights);
    }
};

/*---------------------------------------------------------------------------*\
                        Floor interpolation Declaration
\*---------------------------------------------------------------------------*/

class floor
:
    public interpolationWeight1D
{
public:

    floor(const List<scalar>& xs)
    :
        interpolationWeight1D(xs)
    {}

    virtual ~floor()
    {}

    virtual word type() const
    {
        return "floor";
    }

    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const
    {
        return autoPtr<interpolationWeight1D>(new floor(xs));
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
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

    ceil(const List<scalar>& xs)
    :
        interpolationWeight1D(xs)
    {}

    virtual ~ceil()
    {}

    virtual word type() const
    {
        return "ceil";
    }

    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const
    {
        return autoPtr<interpolationWeight1D>(new ceil(xs));
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
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

    linearExtrapolated(const List<scalar>& xs)
    :
        interpolationWeight1D(xs)
    {}

    virtual ~linearExtrapolated()
    {}

    virtual word type() const
    {
        return "linearExtrapolated";
    }

    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const
    {
        return autoPtr<interpolationWeight1D>(new linearExtrapolated(xs));
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        List<label>& indices,
        List<scalar>& weights
    ) const;
};


class linearClamp
:
    public Clamp<linearExtrapolated>
{
public:

    linearClamp(const List<scalar>& xs)
    :
        Clamp<linearExtrapolated>(xs)
    {}

    virtual ~linearClamp()
    {}

    virtual word type() const
    {
        return "linearClamp";
    }

    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const
    {
        return autoPtr<interpolationWeight1D>(new linearClamp(xs));
    }
};

/*---------------------------------------------------------------------------*\
                        Quadratic interpolation Declaration
\*---------------------------------------------------------------------------*/


class quadraticExtrapolated
:
    public interpolationWeight1D
{
public:

    quadraticExtrapolated(const List<scalar>& xs)
    :
        interpolationWeight1D(xs)
    {}

    virtual ~quadraticExtrapolated()
    {}

    virtual word type() const
    {
        return "quadraticExtrapolated";
    }

    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const
    {
        return autoPtr<interpolationWeight1D>(new quadraticExtrapolated(xs));
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        List<label>& indices,
        List<scalar>& weights
    ) const;
};


class quadraticClamp
:
    public Clamp<quadraticExtrapolated>
{
public:

    quadraticClamp(const List<scalar>& xs)
    :
        Clamp<quadraticExtrapolated>(xs)
    {}

    virtual ~quadraticClamp()
    {}

    virtual word type() const
    {
        return "quadraticClamp";
    }

    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const
    {
        return autoPtr<interpolationWeight1D>(new quadraticClamp(xs));
    }
};


/*---------------------------------------------------------------------------*\
                        Cubic interpolation Declaration
\*---------------------------------------------------------------------------*/


class cubicExtrapolated
:
    public interpolationWeight1D
{
public:

    cubicExtrapolated(const List<scalar>& xs)
    :
        interpolationWeight1D(xs)
    {}

    virtual ~cubicExtrapolated()
    {}

    virtual word type() const
    {
        return "cubicExtrapolated";
    }

    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const
    {
        return autoPtr<interpolationWeight1D>(new cubicExtrapolated(xs));
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        List<label>& indices,
        List<scalar>& weights
    ) const;
};


class cubicClamp
:
    public Clamp<cubicExtrapolated>
{
public:

    cubicClamp(const List<scalar>& xs)
    :
        Clamp<cubicExtrapolated>(xs)
    {}

    virtual ~cubicClamp()
    {}

    virtual word type() const
    {
        return "cubicClamp";
    }

    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const
    {
        return autoPtr<interpolationWeight1D>(new cubicClamp(xs));
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interpolationWeights1D


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
