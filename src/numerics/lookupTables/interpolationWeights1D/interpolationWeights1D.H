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
#include "runTimeSelectionTables.H"
#include "DynamicList.H"
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

    TypeName("interpolationWeight1D");

    interpolationWeight1D(const List<scalar>& xs)
    :
        xs_(xs)
    {}
    interpolationWeight1D(const dictionary& dict, const List<scalar>& xs)
    :
        xs_(xs)
    {}

    declareRunTimeSelectionTable
    (
        autoPtr,
        interpolationWeight1D,
        null,
        (
            const List<scalar>& xs
        ),
        (xs)
    );
    declareRunTimeSelectionTable
    (
        autoPtr,
        interpolationWeight1D,
        dictionary,
        (
            const dictionary& dict,
            const List<scalar>& xs
        ),
        (dict, xs)
    );

    virtual ~interpolationWeight1D()
    {}

    static autoPtr<interpolationWeight1D> New
    (
        const word&,
        const List<scalar>&,
        const bool finalData = true
    );
    static autoPtr<interpolationWeight1D> New
    (
        const word&,
        const dictionary& dict,
        const List<scalar>&,
        const bool finalData = true
    );
    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const = 0;

    virtual bool validate(const bool fail = true) const;
    virtual label nRequired() const = 0;
    virtual void update()
    {}

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
        DynamicList<label>& indices,
        DynamicList<scalar>& weights
    ) const = 0;

    virtual void updateDWeights
    (
        const scalar x,
        const scalar xMod,
        const label i,
        const UList<scalar>& xs,
        DynamicList<label>& indices,
        DynamicList<scalar>& weights,
        DynamicList<scalar>& dws
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
    Clamp(const dictionary& dict, const List<scalar>& xs)
    :
        Interp(dict, xs)
    {}

    virtual ~Clamp()
    {}

    virtual label nRequired() const
    {
        return Interp::nRequired();
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        DynamicList<label>& indices,
        DynamicList<scalar>& weights
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

    virtual void updateDWeights
    (
        const scalar x,
        const scalar xMod,
        const label i,
        const UList<scalar>& xs,
        DynamicList<label>& indices,
        DynamicList<scalar>& weights,
        DynamicList<scalar>& dws
    ) const
    {
        if (x < xs[0])
        {
            indices.setSize(1);
            weights.setSize(1);
            dws.setSize(1);

            indices[0] = 0;
            weights[0] = 0;
            dws[0] = 0.0;
            return;
        }
        else if (x > xs.last())
        {
            indices.setSize(1);
            weights.setSize(1);
            dws.setSize(1);

            indices[0] = xs.size()-1;
            weights[0] = 1.0;
            dws[0] = 0.0;
            return;
        }
        Interp::updateDWeights(x, xMod, i, xs, indices, weights, dws);
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

    TypeName("floor");

    floor(const List<scalar>& xs)
    :
        interpolationWeight1D(xs)
    {}
    floor(const dictionary& dict, const List<scalar>& xs)
    :
        interpolationWeight1D(dict, xs)
    {}

    virtual ~floor()
    {}

    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const
    {
        return autoPtr<interpolationWeight1D>(new floor(xs));
    }

    virtual label nRequired() const
    {
        return 1;
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        DynamicList<label>& indices,
        DynamicList<scalar>& weights
    ) const;

    virtual void updateDWeights
    (
        const scalar x,
        const scalar xMod,
        const label i,
        const UList<scalar>& xs,
        DynamicList<label>& indices,
        DynamicList<scalar>& weights,
        DynamicList<scalar>& dws
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

    TypeName("ceil");

    ceil(const List<scalar>& xs)
    :
        interpolationWeight1D(xs)
    {}
    ceil(const dictionary& dict, const List<scalar>& xs)
    :
        interpolationWeight1D(dict, xs)
    {}

    virtual ~ceil()
    {}

    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const
    {
        return autoPtr<interpolationWeight1D>(new ceil(xs));
    }

    virtual label nRequired() const
    {
        return 1;
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        DynamicList<label>& indices,
        DynamicList<scalar>& weights
    ) const;

    virtual void updateDWeights
    (
        const scalar x,
        const scalar xMod,
        const label i,
        const UList<scalar>& xs,
        DynamicList<label>& indices,
        DynamicList<scalar>& weights,
        DynamicList<scalar>& dws
    ) const;
};


/*---------------------------------------------------------------------------*\
                linearExtrapolated interpolation Declaration
\*---------------------------------------------------------------------------*/\

class linearExtrapolated
:
    public interpolationWeight1D
{
public:

    TypeName("linearExtrapolated");

    linearExtrapolated(const List<scalar>& xs)
    :
        interpolationWeight1D(xs)
    {}
    linearExtrapolated(const dictionary& dict, const List<scalar>& xs)
    :
        interpolationWeight1D(dict, xs)
    {}

    virtual ~linearExtrapolated()
    {}

    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const
    {
        return autoPtr<interpolationWeight1D>(new linearExtrapolated(xs));
    }

    virtual label nRequired() const
    {
        return 2;
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        DynamicList<label>& indices,
        DynamicList<scalar>& weights
    ) const;

    virtual void updateDWeights
    (
        const scalar x,
        const scalar xMod,
        const label i,
        const UList<scalar>& xs,
        DynamicList<label>& indices,
        DynamicList<scalar>& weights,
        DynamicList<scalar>& dws
    ) const;
};


class linearClamp
:
    public Clamp<linearExtrapolated>
{
    bool singular_;

public:

    TypeName("linearClamp");

    linearClamp(const List<scalar>& xs)
    :
        Clamp<linearExtrapolated>(xs),
        singular_(xs.size() < 2)
    {}
    linearClamp(const dictionary& dict, const List<scalar>& xs)
    :
        Clamp<linearExtrapolated>(dict, xs),
        singular_(xs.size() < 2)
    {}

    virtual ~linearClamp()
    {}

    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const
    {
        return autoPtr<interpolationWeight1D>(new linearClamp(xs));
    }

    // This is a special case
    virtual label nRequired() const
    {
        return 1;
    }

    virtual void update()
    {
        singular_ = xs_.size() < 2;
    }

    virtual void updateWeights
    (
        const scalar x,
        const scalar xMod,
        const label i,
        DynamicList<label>& indices,
        DynamicList<scalar>& weights
    ) const
    {
        if (singular_)
        {
            indices(0) = 0;
            weights(0) = 1.0;
        }
        else
        {
            Clamp<linearExtrapolated>::updateWeights(x, i, indices, weights);
        }
    }

    virtual void updateDWeights
    (
        const scalar x,
        const scalar xMod,
        const label i,
        const UList<scalar>& xs,
        DynamicList<label>& indices,
        DynamicList<scalar>& weights,
        DynamicList<scalar>& dws
    ) const
    {
        if (singular_)
        {
            indices(0) = 0;
            weights(0) = 1.0;
            dws(0) = 0.0;
        }
        else
        {
            Clamp<linearExtrapolated>::updateDWeights
            (
                x,
                xMod,
                i,
                xs,
                indices,
                weights,
                dws
            );
        }
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

    TypeName("quadraticExtrapolated");

    quadraticExtrapolated(const List<scalar>& xs)
    :
        interpolationWeight1D(xs)
    {}
    quadraticExtrapolated(const dictionary& dict, const List<scalar>& xs)
    :
        interpolationWeight1D(dict, xs)
    {}

    virtual ~quadraticExtrapolated()
    {}

    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const
    {
        return autoPtr<interpolationWeight1D>(new quadraticExtrapolated(xs));
    }

    virtual label nRequired() const
    {
        return 3;
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        DynamicList<label>& indices,
        DynamicList<scalar>& weights
    ) const;

    virtual void updateDWeights
    (
        const scalar x,
        const scalar xMod,
        const label i,
        const UList<scalar>& xs,
        DynamicList<label>& indices,
        DynamicList<scalar>& weights,
        DynamicList<scalar>& dws
    ) const;
};


class quadraticClamp
:
    public Clamp<quadraticExtrapolated>
{
public:

    TypeName("quadraticClamp");

    quadraticClamp(const List<scalar>& xs)
    :
        Clamp<quadraticExtrapolated>(xs)
    {}
    quadraticClamp(const dictionary& dict, const List<scalar>& xs)
    :
        Clamp<quadraticExtrapolated>(dict, xs)
    {}

    virtual ~quadraticClamp()
    {}

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

    TypeName("cubicExtrapolated");

    cubicExtrapolated(const List<scalar>& xs)
    :
        interpolationWeight1D(xs)
    {}
    cubicExtrapolated(const dictionary& dict, const List<scalar>& xs)
    :
        interpolationWeight1D(dict, xs)
    {}

    virtual ~cubicExtrapolated()
    {}

    virtual autoPtr<interpolationWeight1D> clone(const List<scalar>& xs) const
    {
        return autoPtr<interpolationWeight1D>(new cubicExtrapolated(xs));
    }

    virtual label nRequired() const
    {
        return 4;
    }

    virtual void updateWeights
    (
        const scalar x,
        const label i,
        DynamicList<label>& indices,
        DynamicList<scalar>& weights
    ) const;

    virtual void updateDWeights
    (
        const scalar x,
        const scalar xMod,
        const label i,
        const UList<scalar>& xs,
        DynamicList<label>& indices,
        DynamicList<scalar>& weights,
        DynamicList<scalar>& dws
    ) const;
};


class cubicClamp
:
    public Clamp<cubicExtrapolated>
{
public:

    TypeName("cubicClamp");

    cubicClamp(const List<scalar>& xs)
    :
        Clamp<cubicExtrapolated>(xs)
    {}
    cubicClamp(const dictionary& dict, const List<scalar>& xs)
    :
        Clamp<cubicExtrapolated>(dict, xs)
    {}

    virtual ~cubicClamp()
    {}

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
