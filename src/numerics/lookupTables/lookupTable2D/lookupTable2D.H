/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2022
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
    Foam::lookupTable2D

Description
    Table used to lookup values given a 2D table

SourceFiles
    lookupTable2D.C

\*---------------------------------------------------------------------------*/

#ifndef lookupTable2D_H
#define lookupTable2D_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "List2D.H"
#include "IOstreams.H"
#include "Switch.H"
#include "Field.H"
#include "fileName.H"

#include "Modifier.H"
#include "interpolationWeights1D.H"
#include "indexing.H"

#include "EquationsFwd.H"
#include "UnivariateEquationsFwd.H"
#include "univariateRootSolver.H"
#include "NewtonRaphsonUnivariateRootSolver.H"

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Type> class lookupTable2D;
template<class Type> Ostream& operator<<(Ostream&, const lookupTable2D<Type>&);

/*---------------------------------------------------------------------------*\
                           Class lookupTable2D Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class lookupTable2D
{
public:

    typedef Vector2D<label> labelVector2D;

protected:
// Protected data

    //- Names of components
    word xName_;
    word yName_;
    word fName_;

    //- Modifiers
    autoPtr<Modifier<Type>> mod_;
    autoPtr<Modifier<scalar>> modX_;
    autoPtr<Modifier<scalar>> modY_;

    //- Data
    List2D<Type> data_;

    //- Modified x field values
    List<scalar> xModValues_;

    //- Modified y field values
    List<scalar> yModValues_;

    //- Lookup indexes in the x direction
    autoPtr<indexer> xIndexing_;

    //- Lookup indexes in the y direction
    autoPtr<indexer> yIndexing_;

    //- X-interpolater type
    autoPtr<interpolationWeight1D> xInterpolator_;

    //- Y-interpolater type
    autoPtr<interpolationWeight1D> yInterpolator_;

    //- Stored real data
    List2D<Type>* realDataPtr_;

    //- Stored real x values
    List<scalar>* xValuesPtr_;

    //- Stored real y values
    List<scalar>* yValuesPtr_;

    //- Current indexes (used to access location in list after lookup)
    mutable labelVector2D ij_;
    mutable DynamicList<label> is_;
    mutable DynamicList<label> js_;
    mutable DynamicList<labelVector2D> indices_;

    //- Current interpolation weights
    mutable DynamicList<scalar> weights_;
    mutable DynamicList<scalar> wxs_;
    mutable DynamicList<scalar> wys_;

    // Equation to evaulate the table
    class equation2D
    :
        public UnivariateEquation<Type>
    {
        const lookupTable2D<Type>& table_;

    public:

        equation2D(const lookupTable2D<Type>& table)
        :
            UnivariateEquation<Type>(2),
            table_(table)
        {}

        label nDerivatives() const
        {
            return 1;
        }

        inline Type fx
        (
            const typename UnivariateEquation<Type>::XType& x,
            const label li
        ) const
        {
            return table_.lookup(x[0], x[1]);
        }

        inline void dfdX
        (
            const typename UnivariateEquation<Type>::XType& x,
            const label li,
            List<Type>& dfdx
        ) const
        {
            dfdx[0] = table_.dFdX(x[0], x[1]);
            dfdx[1] = table_.dFdY(x[0], x[1]);
        }
    };

    // Equation to used to find a give component
    class errorEquation
    :
        public Equation<Type>
    {
        const lookupTable2D<Type>& table_;
        label dir_;
        scalar y_;
        Type f_;

    public:

        errorEquation(const lookupTable2D<Type>& table)
        :
            Equation<Type>(),
            table_(table),
            dir_(-1)
        {}

        void setDir(const label dir)
        {
            if (dir == dir_)
            {
                return;
            }
            dir_ = dir;
            if (dir == 0)
            {
                this->setLower(min(table_.xValues()));
                this->setUpper(max(table_.xValues()));
            }
            else if (dir == 1)
            {
                this->setLower(min(table_.yValues()));
                this->setUpper(max(table_.yValues()));
            }
            else
            {
                FatalErrorInFunction
                    << "Table has 2 dimensions, but trying to find dimension "
                    << dir+1 << endl
                    << abort(FatalError);
            }
        }

        void set(const label dir, const Type& f, const scalar y)
        {
            y_ = y;
            f_ = f;
            setDir(dir);
        }
        label nDerivatives() const
        {
            return 2;
        }

        virtual Type fx(const scalar x, const label li) const
        {
            if (dir_ == 0)
            {
                return table_.lookup(x, y_) - f_;
            }
            else if (dir_ == 1)
            {
                return table_.lookup(y_, x) - f_;
            }
            else
            {
                FatalErrorInFunction
                    << "direction has not been set" << endl
                    << abort(FatalError);
            }
            return Zero;
        }

        virtual Type dfdx(const scalar x, const label li) const
        {
            if (dir_ == 0)
            {
                return table_.dFdX(x, y_);
            }
            else if (dir_ == 1)
            {
                return table_.dFdY(y_, x);
            }
            else
            {
                FatalErrorInFunction
                    << "direction has not been set" << endl
                    << abort(FatalError);
            }
            return Zero;
        }

        virtual Type d2fdx2(const scalar x, const label li) const
        {
            if (dir_ == 0)
            {
                return table_.d2FdX2(x, y_);
            }
            else if (dir_ == 1)
            {
                return table_.d2FdY2(y_, x);
            }
            else
            {
                FatalErrorInFunction
                    << "direction has not been set" << endl
                    << abort(FatalError);
            }
            return Zero;
        }
    };


    //- Error equation
    mutable autoPtr<errorEquation> error_;

    //- Root solver
    mutable autoPtr<univariateRootSolver> solver_;

    // Protected member functions

        //- Return the real data
        inline const List2D<Type>& fValues() const
        {
            return *realDataPtr_;
        }

        //- Access the real data values
        inline List2D<Type>& fValues()
        {
            return *realDataPtr_;
        }

        //- Return the real x values
        inline const List<scalar>& xValues() const
        {
            return *xValuesPtr_;
        }

        //- Access the real x values
        inline List<scalar>& xValues()
        {
            return *xValuesPtr_;
        }

        //- Return the real y values
        inline const List<scalar>& yValues() const
        {
            return *yValuesPtr_;
        }

        //- Access the real y values
        inline List<scalar>& yValues()
        {
            return *yValuesPtr_;
        }


public:

    // Constructors

        //- Construct empty
        lookupTable2D();

        //- Copy constructor
        lookupTable2D(const lookupTable2D&);

        //- Construct from dictionary
        lookupTable2D
        (
            const dictionary& dict,
            const word& xName,
            const word& yName,
            const word& name,
            const bool canRead = true
        );

        //- Construct from data
        template<template<class> class ListListType>
        lookupTable2D
        (
            const List<scalar>& x,
            const List<scalar>& y,
            const ListListType<Type>& data,
            const word& xMod,
            const word& yMod,
            const word& mod,
            const word& xInterpolationScheme,
            const word& yInterpolationScheme,
            const bool isReal = true
        );



    //- Destructor
    virtual ~lookupTable2D();

    // Member functions

        // Access

            //- Return x name
            const word& xName() const
            {
                return xName_;
            }

            //- Access x name
            word& xName()
            {
                return xName_;
            }

            //- Return y name
            const word& yName() const
            {
                return yName_;
            }

            //- Access y name
            word& yName()
            {
                return yName_;
            }

            //- Return f name
            const word& fName() const
            {
                return fName_;
            }

            //- Access f name
            word& fName()
            {
                return fName_;
            }

            //- Return current index (x)
            label i() const
            {
                return ij_.x();
            }

            //- Return current index (y)
            label j() const
            {
                return ij_.y();
            }

            //- Return current interpolation indices
            const List<labelVector2D>& indices() const
            {
                return indices_;
            }

            //- Return current interpolation indices
            const scalarList& weights() const
            {
                return weights_;
            }

            //- Return the modifier
            const Modifier<Type>& mod() const
            {
                return mod_();
            }

            //- Return the x-modifier
            const Modifier<scalar>& modX() const
            {
                return modX_();
            }

            //- Return the y-modifier
            const Modifier<scalar>& modY() const
            {
                return modY_();
            }

            //- Return to real x values
            const List<scalar>& x() const
            {
                return xValues();
            }

            //- Return to real y values
            const List<scalar>& y() const
            {
                return yValues();
            }

            //- Return to modified x values
            const List<scalar>& xMod() const
            {
                return xModValues_;
            }

            //- Return to modified y values
            const List<scalar>& yMod() const
            {
                return yModValues_;
            }

            //- Return to modified data values
            const List2D<Type>& fMod() const
            {
                return data_;
            }

            //- Return real data values
            const List2D<Type>& f() const
            {
                return fValues();
            }

            //- Return the interpolation scheme for x
            const interpolationWeight1D& xInterpolator() const
            {
                return xInterpolator_();
            }

            //- Return the interpolation scheme for y
            const interpolationWeight1D& yInterpolator() const
            {
                return yInterpolator_();
            }


        // Edit functions

            //- Update indexes
            void updateIndex(const scalar x, const scalar y) const;

            //- Update indexes and weights
            void update(const scalar x, const scalar y) const;

            //- Lookup value
            Type lookup(const scalar x, const scalar y) const;

            //- Lookup x given f and y
            scalar reverseLookupX(const Type& fin, const scalar y) const;

            //- Lookup y given f and x
            scalar reverseLookupY(const Type& fin, const scalar x) const;

            //- Interpolate value using the current index and weight
            template<template<class> class ListListType, class fType>
            fType interpolate(const ListListType<fType>&) const;

            //- Interpolate value using the current index and weight
            template<template<class> class ListListType, class fType>
            fType interpolate
            (
                const scalar,
                const scalar,
                const ListListType<fType>&
            ) const;

            //- Return first derivative w.r.t. x
            Type dFdX(const scalar x, const scalar y) const;

            //- Return first derivative w.r.t y
            Type dFdY(const scalar x, const scalar y) const;

            //- Return second derivative w.r.t. x
            Type d2FdX2(const scalar x, const scalar y) const;

            //- Return mixed second derivative
            Type d2FdXdY(const scalar x, const scalar y) const;

            //- Return second derivative w.r.t. y
            Type d2FdY2(const scalar x, const scalar y) const;

            //- Return value from a given list from lower index and weight
            Type getValue(const label ij, const scalar f, const UList<Type>& xy) const;


            //- Set data, x, and y
            template<template<class> class ListListType>
            void set
            (
                const List<scalar>& x,
                const List<scalar>& y,
                const ListListType<Type>& data,
                const bool isReal = true
            );
            template<template<class> class ListListType>
            void set
            (
                const List<scalar>& x,
                const List<scalar>& y,
                const ListListType<Type>& data,
                const word& xMod,
                const word& yMod,
                const word& mod,
                const word& xInterpolationScheme,
                const word& yInterpolationScheme,
                const bool isReal = true
            );

            //- Set x
            void setX
            (
                const List<scalar>& x,
                const bool isReal = true
            );
            void setX
            (
                const List<scalar>& x,
                const word& xMod,
                const bool isReal = true
            );

            //- Set y
            void setY
            (
                const List<scalar>& y,
                const bool isReal = true
            );
            void setY
            (
                const List<scalar>& y,
                const word& yMod,
                const bool isReal = true
            );

            //- Set data
            template<template<class> class ListListType>
            void setData
            (
                const ListListType<Type>& data,
                const bool isReal = true
            );
            template<template<class> class ListListType>
            void setData
            (
                const ListListType<Type>& data,
                const word& mod,
                const bool isReal = true
            );

        //- Return an equation using the table
        autoPtr<univariateEquation<Type>> eqn() const
        {
            return autoPtr<univariateEquation<Type>>(new equation2D(*this));
        }

        //- Return an equation using the table
        errorEquation& error(const label dir = -1) const
        {
            if (!error_.valid())
            {
                error_.set(new errorEquation(*this));
            }
            if (dir != -1)
            {
                error_->setDir(dir);
            }
            return error_();
        }

        //- Return an equation using the table
        const univariateRootSolver& solver
        (
            const word& type = NewtonRaphsonUnivariateRootSolver::typeName,
            const dictionary& dict = dictionary()
        ) const
        {
            if (!solver_.valid())
            {
                solver_ = univariateRootSolver::New
                (
                    type,
                    error(),
                    dict
                );
            }
            return solver_();
        }
        //- Return an equation using the table
        const univariateRootSolver& solver
        (
            const label dir,
            const Type& f,
            const scalar y,
            const word& type = NewtonRaphsonUnivariateRootSolver::typeName,
            const dictionary& dict = dictionary()
        ) const
        {
            error().set(dir, f, y);
            return solver(type, dict);
        }

        //- Read from a given dictionary
        void read
        (
            const dictionary& dict,
            const word& xName,
            const word& yName,
            const word& name,
            const bool canRead = true
        );

    //- Write data to dictionary stream
    virtual void write(Ostream& os, const word& dictName = word::null) const;


    // Member Operators

        //- Assignment
        void operator=(const lookupTable2D<Type>&);


    // IOstream Operators

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const lookupTable2D<Type>& table
        );
};


template<class Type>
void writeEntry(Ostream& os, const lookupTable2D<Type>& table);

template<class Type>
void writeEntry
(
    Ostream& os,
    const word& dictName,
    const lookupTable2D<Type>& table
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

template<>
Foam::scalar Foam::lookupTable2D<Foam::scalar>::reverseLookupX
(
    const scalar& fin,
    const scalar y
) const;

template<>
Foam::scalar Foam::lookupTable2D<Foam::scalar>::reverseLookupY
(
    const scalar& fin,
    const scalar x
) const;

template<>
void Foam::lookupTable2D<Foam::scalar>::read
(
    const dictionary& dict,
    const word& xName,
    const word& yName,
    const word& name,
    const bool canRead
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lookupTable2D.C"
    #include "lookupTable2DTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
