/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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
    Foam::lookupTable3D

Description
    Table used to lookup values given a 3D table

SourceFiles
    lookupTable3D.C

\*---------------------------------------------------------------------------*/

#ifndef lookupTable3D_H
#define lookupTable3D_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "List3D.H"
#include "IOstreams.H"
#include "Switch.H"
#include "fileName.H"
#include "labelVector.H"

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
template<class Type> class lookupTable3D;
template<class Type> Ostream& operator<<(Ostream&, const lookupTable3D<Type>&);


/*---------------------------------------------------------------------------*\
                           Class lookupTable3D Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class lookupTable3D
{
protected:
// Protected data

    //- Names of components
    word xName_;
    word yName_;
    word zName_;
    word fName_;

    //- Modifiers
    autoPtr<Modifier<Type>> mod_;
    autoPtr<Modifier<scalar>> modX_;
    autoPtr<Modifier<scalar>> modY_;
    autoPtr<Modifier<scalar>> modZ_;

    //- Data
    List3D<Type> data_;

    //- Modified x field values
    List<scalar> xModValues_;

    //- Modified y field values
    List<scalar> yModValues_;

    //- Modified z field values
    List<scalar> zModValues_;

    //- Lookup indexes in the x direction
    autoPtr<indexer> xIndexing_;

    //- Lookup indexes in the y direction
    autoPtr<indexer> yIndexing_;

    //- Lookup indexes in the z direction
    autoPtr<indexer> zIndexing_;

    //- X-interpolator
    autoPtr<interpolationWeight1D> xInterpolator_;

    //- Y-interpolator
    autoPtr<interpolationWeight1D> yInterpolator_;

    //- Z-interpolator
    autoPtr<interpolationWeight1D> zInterpolator_;

    //- Stored real data
    List3D<Type>* realDataPtr_;

    //- Stored real x values
    List<scalar>* xValuesPtr_;

    //- Stored real y values
    List<scalar>* yValuesPtr_;

    //- Stored real z values
    List<scalar>* zValuesPtr_;

    //- Current indexes (used to access location in list after lookup)
    mutable labelVector ijk_;
    mutable DynamicList<label> is_;
    mutable DynamicList<label> js_;
    mutable DynamicList<label> ks_;
    mutable DynamicList<labelVector> indices_;

    //- Current interpolation weights
    mutable DynamicList<scalar> wxs_;
    mutable DynamicList<scalar> wys_;
    mutable DynamicList<scalar> wzs_;
    mutable DynamicList<scalar> weights_;


    // Equation to evaulate the table
    class equation3D
    :
        public UnivariateEquation<Type>
    {
        const lookupTable3D<Type>& table_;

    public:

        equation3D(const lookupTable3D<Type>& table)
        :
            UnivariateEquation<Type>(3),
            table_(table)
        {
            this->setLowerLimits
            (
                {
                    min(table_.xValues()),
                    min(table_.yValues()),
                    min(table_.zValues())
                }
            );
            this->setUpperLimits
            (
                {
                    max(table_.xValues()),
                    max(table_.yValues()),
                    max(table_.zValues())
                }
            );
        }

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
            return table_.lookup(x[0], x[1], x[2]);
        }

        inline void dfdX
        (
            const typename UnivariateEquation<Type>::XType& x,
            const label li,
            List<Type>& dfdx
        ) const
        {
            dfdx[0] = table_.dFdX(x[0], x[1], x[2]);
            dfdx[1] = table_.dFdY(x[0], x[1], x[2]);
            dfdx[2] = table_.dFdZ(x[0], x[1], x[2]);
        }
    };

    // Equation to used to find a give component
    class errorEquation
    :
        public Equation<Type>
    {
        const lookupTable3D<Type>& table_;
        label dir_;
        scalar y_;
        scalar z_;
        Type f_;

    public:

        errorEquation(const lookupTable3D<Type>& table)
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
            else if (dir == 2)
            {
                this->setLower(min(table_.zValues()));
                this->setUpper(max(table_.zValues()));
            }
            else
            {
                FatalErrorInFunction
                    << "Table has 3 dimensions, but trying to find dimension "
                    << dir+1 << endl
                    << abort(FatalError);
            }
        }

        void set(const label dir, const Type& f, const scalar y, const scalar z)
        {
            y_ = y;
            z_ = z;
            f_ = f;
            setDir(dir);
        }
        label nDerivatives() const
        {
            return 2;
        }

        virtual scalar fx(const scalar x, const label li) const
        {
            if (dir_ == 0)
            {
                return table_.lookup(x, y_, z_) - f_;
            }
            else if (dir_ == 1)
            {
                return table_.lookup(y_, x, z_) - f_;
            }
            else if (dir_ == 2)
            {
                return table_.lookup(y_, z_, x) - f_;
            }
            else
            {
                FatalErrorInFunction
                    << "direction has not been set" << endl
                    << abort(FatalError);
            }
            return 0.0;
        }

        virtual scalar dfdx(const scalar x, const label li) const
        {
            if (dir_ == 0)
            {
                return table_.dFdX(x, y_, z_);
            }
            else if (dir_ == 1)
            {
                return table_.dFdY(y_, x, z_);
            }
            else if (dir_ == 2)
            {
                return table_.dFdZ(y_, z_, x);
            }
            else
            {
                FatalErrorInFunction
                    << "direction has not been set" << endl
                    << abort(FatalError);
            }
            return 0.0;
        }

        virtual scalar d2fdx2(const scalar x, const label li) const
        {
            if (dir_ == 0)
            {
                return table_.d2FdX2(x, y_, z_);
            }
            else if (dir_ == 1)
            {
                return table_.d2FdY2(y_, x, z_);
            }
            else if (dir_ == 2)
            {
                return table_.d2FdZ2(y_, z_, x);
            }
            else
            {
                FatalErrorInFunction
                    << "direction has not been set" << endl
                    << abort(FatalError);
            }
            return 0.0;
        }
    };


    //- Error equation
    mutable autoPtr<errorEquation> error_;

    //- Root solver
    mutable autoPtr<univariateRootSolver> solver_;

    // Protected member functions

        //- Return the real data
        inline const List3D<scalar>& fValues() const
        {
            return *realDataPtr_;
        }

        //- Access the real data values
        inline List3D<scalar>& fValues()
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

        //- Return the real z values
        inline const List<scalar>& zValues() const
        {
            return *zValuesPtr_;
        }

        //- Access the real z values
        inline List<scalar>& zValues()
        {
            return *zValuesPtr_;
        }


public:

    // Constructors

        //- Construct empty
        lookupTable3D();

        //- Copy constructor
        lookupTable3D(const lookupTable3D&);

        //- Construct from dictionary
        lookupTable3D
        (
            const dictionary& dict,
            const word& xName,
            const word& yName,
            const word& zName,
            const word& name,
            const bool canRead = true
        );

        //- Construct from data
        template<class ListListListType>
        lookupTable3D
        (
            const List<scalar>& x,
            const List<scalar>& y,
            const List<scalar>& z,
            const ListListListType& data,
            const word& xMod,
            const word& yMod,
            const word& zMod,
            const word& mod,
            const word& xInterpolationScheme,
            const word& yInterpolationScheme,
            const word& zInterpolationScheme,
            const bool isReal = true
        );



    //- Destructor
    virtual ~lookupTable3D();


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

            //- Return z name
            const word& zName() const
            {
                return zName_;
            }

            //- Access z name
            word& zName()
            {
                return zName_;
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
                return ijk_.x();
            }

            //- Return current index (y)
            label j() const
            {
                return ijk_.y();
            }

            //- Return current index (z)
            label k() const
            {
                return ijk_.z();
            }

            //- Return current interpolation indices
            const List<labelVector>& indices() const
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

            //- Return the z-modifier
            const Modifier<scalar>& modZ() const
            {
                return modZ_();
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

            //- Return to real z values
            const List<scalar>& z() const
            {
                return zValues();
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

            //- Return to modified z values
            const List<scalar>& zMod() const
            {
                return zModValues_;
            }

            //- Return to modified data values
            const List3D<Type>& fMod() const
            {
                return data_;
            }

            //- Return real data values
            const List3D<Type>& f() const
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

            //- Return the interpolation scheme for z
            const interpolationWeight1D& zInterpolator() const
            {
                return zInterpolator_();
            }


        // Lookup Functions

            //- Update indexes
            void updateIndex(const scalar x, const scalar y, const scalar z) const;

            //- Update indexes and weights
            void update(const scalar x, const scalar y, const scalar z) const;

            //- Lookup value
            Type lookup(const scalar x, const scalar y, const scalar z) const;

            //- Lookup x given f and y and z
            scalar reverseLookupX(const Type& fin, const scalar y, const scalar z) const;

            //- Lookup y given f and z and z
            scalar reverseLookupY(const Type& fin, const scalar x, const scalar z) const;

            //- Lookup z given f and x and y
            scalar reverseLookupZ(const Type& fin, const scalar x, const scalar y) const;

            //- Return first derivative w.r.t.. x
            Type dFdX(const scalar x, const scalar y, const scalar z) const;

            //- Return first derivative w.r.t. y
            Type dFdY(const scalar x, const scalar y, const scalar z) const;

            //- Return first derivative w.r.t. z
            Type dFdZ(const scalar x, const scalar y, const scalar z) const;

            //- Return second derivative w.r.t. x
            Type d2FdX2(const scalar x, const scalar y, const scalar z) const;

            //- Return second derivative w.r.t. y
            Type d2FdY2(const scalar x, const scalar y, const scalar z) const;

            //- Return second derivative w.r.t. z
            Type d2FdZ2(const scalar x, const scalar y, const scalar z) const;

            //- Return mixed second derivative w.r.t. x and y
            Type d2FdXdY(const scalar x, const scalar y, const scalar z) const;

            //- Return mixed second derivative  w.r.t. x and z
            Type d2FdXdZ(const scalar x, const scalar y, const scalar z) const;

            //- Return mixed second derivative  w.r.t. y and z
            Type d2FdYdZ(const scalar x, const scalar y, const scalar z) const;

            //- Return value from a given list from lower index and weight
            Type getValue
            (
                const label ijk,
                const scalar f,
                const List<Type>& xyz
            ) const;

            //- Interpolate value using the current index and weight
            template<class fType, template<class> class ListListList>
            fType interpolate(const ListListList<fType>&) const;

            //- Interpolate value using the current index and weight
            template<class fType, template<class> class ListListList>
            fType interpolate
            (
                const scalar,
                const scalar,
                const scalar,
                const ListListList<fType>&
            ) const;


        // Edit functions

            //- Set data, x, and y
            template<class ListListListType>
            void set
            (
                const List<scalar>& x,
                const List<scalar>& y,
                const List<scalar>& z,
                const ListListListType& data,
                const bool isReal = true
            );
            template<class ListListListType>
            void set
            (
                const List<scalar>& x,
                const List<scalar>& y,
                const List<scalar>& z,
                const ListListListType& data,
                const word& xMod,
                const word& yMod,
                const word& zMod,
                const word& mod,
                const word& xInterpolationScheme,
                const word& yInterpolationScheme,
                const word& zInterpolationScheme,
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

            //- Set z
            void setZ
            (
                const List<scalar>& z,
                const bool isReal = true
            );
            void setZ
            (
                const List<scalar>& z,
                const word& zMod,
                const bool isReal = true
            );

            //- Set data
            template<class ListListListType>
            void setData
            (
                const ListListListType& data,
                const bool isReal = true
            );
            template<class ListListListType>
            void setData
            (
                const ListListListType& data,
                const word& mod,
                const bool isReal = true
            );

        //- Return an equation using the table
        autoPtr<univariateEquation<Type>> eqn() const
        {
            return autoPtr<univariateEquation<Type>>(new equation3D(*this));
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
            const scalar z,
            const word& type = NewtonRaphsonUnivariateRootSolver::typeName,
            const dictionary& dict = dictionary()
        ) const
        {
            error().set(dir, f, y, z);
            return solver(type, dict);
        }


        //- Read from a given dictionary
        void read
        (
            const dictionary& dict,
            const word& xName,
            const word& yName,
            const word& zName,
            const word& name,
            const bool canRead = true
        );

    //- Write data to dictionary stream
    virtual void write(Ostream& os) const;


    // Member Operators

        //- Assignment
        void operator=(const lookupTable3D<Type>&);


    // IOstream Operators

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const lookupTable3D<Type>& table
        );
};


template<class Type>
void writeEntry(Ostream& os, const lookupTable3D<Type>& table);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


template<>
Foam::scalar Foam::lookupTable3D<Foam::scalar>::reverseLookupX
(
    const scalar& fin,
    const scalar y,
    const scalar z
) const;

template<>
Foam::scalar Foam::lookupTable3D<Foam::scalar>::reverseLookupY
(
    const scalar& fin,
    const scalar x,
    const scalar y
) const;

template<>
Foam::scalar Foam::lookupTable3D<Foam::scalar>::reverseLookupZ
(
    const scalar& fin,
    const scalar y,
    const scalar z
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lookupTable3D.C"
    #include "lookupTable3DTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
