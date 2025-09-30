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
    Foam::PtrList2D

Description
    2D pointer List

SourceFiles
    PtrList2D.C

\*---------------------------------------------------------------------------*/

#ifndef PtrList2D_H
#define PtrList2D_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "List2D.H"
#include "Vector2D.H"
#include "autoPtr.H"
#include "tmp.H"

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Type> class PtrList2D;

template<class Type> Ostream& operator<<(Ostream&, const PtrList2D<Type>&);
// template<class Type> Istream& operator>>(Istream&, PtrList2D<Type>&);

/*---------------------------------------------------------------------------*\
                           Class List2D Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PtrList2D
:
    public PtrList<Type>
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
        inline PtrList2D();

        //- Construct from size
        PtrList2D(const label m, const label n);

        //- Construct from data and sizes
        PtrList2D(Type** v, const label m, const label n);

        //- Copy constructor
        PtrList2D(const PtrList2D<Type>& lst);

        //- Move constructor
        PtrList2D(PtrList2D<Type>&& lst);

        //- Istream constructor
        PtrList2D(Istream& is);



    //- Destructor
    ~PtrList2D();


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

            //- Set the sizes
            void setSize(const label m, const label n);

            //- Set element i, j
            inline autoPtr<Type> set(const labelVector2D&, Type* v);

            //- Set element i, j
            inline autoPtr<Type> set(const label i, const label j, Type* v);

            //- Set element i, j
            inline autoPtr<Type> set(const labelVector2D&, autoPtr<Type>& v);

            //- Set element i, j
            inline autoPtr<Type> set(const label i, const label j, autoPtr<Type>& v);

            //- Set element i, j
            inline autoPtr<Type> set(const labelVector2D&, tmp<Type>& v);

            //- Set element i, j
            inline autoPtr<Type> set(const label i, const label j, tmp<Type>& v);


            //- Transfer to another PtrList2D
            void transfer(PtrList2D<Type>& lst);


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
        inline UPtrList<Type> operator()(const label) const;

        //- Access element row i
        inline UPtrList<Type> operator()(const label);

        //- Return element row i
        inline UPtrList<Type> operator[](const label) const;

        //- Access element row i
        inline UPtrList<Type> operator[](const label);

        //- Assignment
        void operator=(const PtrList2D<Type>&);


    // IOstream Operators

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const PtrList2D<Type>& table
        );

        // //- Istream Operator
        // friend Istream& operator>> <Type>
        // (
        //     Istream& is,
        //     PtrList2D<Type>& table
        // );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PtrList2DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PtrList2D.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
