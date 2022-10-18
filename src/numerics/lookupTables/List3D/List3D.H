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
    Foam::List3D

Description
    3D List

SourceFiles
    List3D.C

\*---------------------------------------------------------------------------*/

#ifndef List3D_H
#define List3D_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "List.H"
#include "List2D.H"
#include "Vector.H"

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Type> class List3D;

template<class Type> Ostream& operator<<(Ostream&, const List3D<Type>&);
template<class Type> Istream& operator>>(Istream&, List3D<Type>&);

/*---------------------------------------------------------------------------*\
                           Class List3D Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class List3D
:
    public List<Type>
{
public:
    // Public typedef
    typedef Vector<label> labelVector;


protected:
// Protected data

    //- Size in the first direction
    label m_;

    //- Size in the second direction
    label n_;

    //- Size in the third direction
    label l_;

    //- Is the data owned by this object
    bool own_;


    // Proctected memeber functions

        //- Check if the index is valid
        inline bool checkIIndex(const label i) const;

        //- Check if the second index is valid
        inline bool checkJIndex(const label j) const;

        //- Check if the third index is valid
        inline bool checkKIndex(const label k) const;


public:

    // Constructors

        //- Null construct
        inline List3D();

        //- Construct from size
        List3D(const label m, const label n, const label l);

        //- Construct from size and intial value
        List3D(const label m, const label n, const label l, const Type& t);

        //- Construct from size and intial value
        template<class Type2>
        List3D(const label m, const label n, const label l, const Type2& t);

        //- Construct from data and sizes
        List3D(Type* v, const label m, const label n, const label l);

        //- Copy constructor
        List3D(const List3D& lst);

        //- Move constructor
        List3D(List3D<Type>&& lst);

        //- Istream constructor
        List3D(Istream& is);



    //- Destructor
    ~List3D();


    // Member functions

        // Access

            //- Return the 1D index
            inline label index(const label, const label, const label) const;

            //- Return the size in the first direction
            inline label m() const;

            //- Return the size in the second direction
            inline label n() const;

            //- Return the size in the third direction
            inline label l() const;

            //- Return the full size
            inline label size() const;


        // Edit

            //- Set the sizes without updating current entries
            void resize(const label m, const label n, const label l);

            //- Set the sizes and initialize without updating current entries
            void resize
            (
                const label m,
                const label n,
                const label l,
                const Type& t
            );

            //- Set the sizes
            void setSize(const label m, const label n, const label l);

            //- Set the sizes and initialize new values
            void setSize
            (
                const label m,
                const label n,
                const label l,
                const Type& t
            );

            //- Transfer to another List3D
            void transfer(List3D<Type>& lst);


    // Member Operators

        //- Return element i, j, k
        inline const Type& operator()(const labelVector&) const;

        //- Return element i, j, k
        inline Type& operator()(const labelVector&);

        //- Return element i, j, k
        inline const Type& operator()
        (
            const label,
            const label,
            const label
        ) const;

        //- Access element i, j, k
        inline Type& operator()(const label, const label, const label);

        //- Return row k
        inline UList<Type> operator()(const label, const label) const;

        //- Access row k
        inline UList<Type> operator()(const label, const label);

        //- Return element row i
        inline List2D<Type> operator()(const label) const;

        //- Access element row i
        inline List2D<Type> operator()(const label);

        //- Return element row i
        inline List2D<Type> operator[](const label) const;

        //- Access element row i
        inline List2D<Type> operator[](const label);

        //- Assignment
        void operator=(const List3D<Type>&);

        //- Assignment
        template<class Type2>
        void operator=(const List3D<Type2>&);

        //- Assignment
        template<class ListListListType>
        void operator=(const ListListListType&);


    // IOstream Operators

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const List3D<Type>& table
        );

        //- Istream Operator
        friend Istream& operator>> <Type>
        (
            Istream& is,
            List3D<Type>& table
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "List3DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "List3D.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
