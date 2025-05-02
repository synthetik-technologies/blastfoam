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
    Foam::List2D

Description
    2D List

SourceFiles
    List2D.C

\*---------------------------------------------------------------------------*/

#ifndef List2D_H
#define List2D_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "List.H"
#include "Vector2D.H"

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Type> class List2D;

template<class Type> Ostream& operator<<(Ostream&, const List2D<Type>&);
template<class Type> Istream& operator>>(Istream&, List2D<Type>&);

/*---------------------------------------------------------------------------*\
                           Class List2D Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class List2D
:
    public List<Type>
{
public:
    // Public typedef
    typedef Vector2D<label> labelVector2D;


protected:
// Protected data

    //- Size in the first direction
    label m_;

    //- Size in the second direction
    label n_;

    //- Is the data owned by this object
    bool own_;

    // Proctected memeber functions

        //- Check if the index is valid
        inline bool checkIIndex(const label i) const;

        //- Check if the second index is valid
        inline bool checkJIndex(const label j) const;


public:

    // Constructors

        //- Null construct
        inline List2D();

        //- Construct from size
        List2D(const label m, const label n);

        //- Construct from size and intial value
        List2D(const label m, const label n, const Type& t);

        //- Construct from size and intial value
        template<class Type2>
        List2D(const label m, const label n, const Type2& t);

        //- Construct from data and sizes
        List2D(Type* v, const label m, const label n);

        //- Copy constructor
        List2D(const List2D<Type>& lst);

        //- Copy constructor
        template<class ListListType>
        List2D(const ListListType& lst);

        //- Move constructor
        List2D(List2D<Type>&& lst);

        //- Istream constructor
        List2D(Istream& is);



    //- Destructor
    ~List2D();


    // Member functions

        // Access

            //- Return the 1D index
            inline label index(const label, const label) const;

            //- Return the size in the first direction
            inline label m() const;

            //- Return the size in the second direction
            inline label n() const;

            //- Return the full size
            inline label size() const;


        // Edit

            //- Set the sizes without updating current entries
            void resize(const label m, const label n);

            //- Set the sizes and initialize without updating current entries
            void resize(const label m, const label n, const Type& t);

            //- Set the sizes
            void setSize(const label m, const label n);

            //- Set the sizes and initialize new values
            void setSize(const label m, const label n, const Type& t);

            //- Transfer to another List2D
            void transfer(List2D<Type>& lst);

            //- Flip table (transpose)
            void flip();



    // Member Operators

        //- Return element i, j
        inline const Type& operator()(const labelVector2D&) const;

        //- Return element i, j
        inline Type& operator()(const labelVector2D&);

        //- Return element i, j
        inline const Type& operator()(const label, const label) const;

        //- Access element i, j
        inline Type& operator()(const label, const label);

        //- Return element row i
        inline UList<Type> operator()(const label) const;

        //- Access element row i
        inline UList<Type> operator()(const label);

        //- Return element row i
        inline UList<Type> operator[](const label) const;

        //- Access element row i
        inline UList<Type> operator[](const label);

        //- Assignment
        void operator=(const List2D<Type>&);

        //- Assignment
        template<class Type2>
        void operator=(const List2D<Type2>&);

        //- Assignment
        template<class ListListType>
        void operator=(const ListListType&);


    // IOstream Operators

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const List2D<Type>& table
        );

        //- Istream Operator
        friend Istream& operator>> <Type>
        (
            Istream& is,
            List2D<Type>& table
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "List2DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "List2D.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
