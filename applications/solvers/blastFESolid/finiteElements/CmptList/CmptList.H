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
    Foam::CmptList

Description
    Class to access a component of a list of type T

SourceFiles
    CmptListI.H

\*---------------------------------------------------------------------------*/

#ifndef CmptList_H
#define CmptList_H

#include "UList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class CmptList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class CmptList
:
    public UList<T>
{
    //- Component index
    const label cmpt_;


public:

    // Constructors

        //- Construct from list non-const reference and component index
        inline CmptList
        (
            UList<T>& lst,
            const label cmpt
        );

        //- Construct from list const reference and component index
        inline CmptList
        (
            const UList<T>& lst,
            const label cmpt
        );


    // Member functions

        // Access

            //- Return the first element of the list
            inline scalar& first();

            //- Return first element of the list
            inline const scalar& first() const;

            //- Return the last element of the list
            inline scalar& last();

            //- Return the last element of the list
            inline const scalar& last() const;


    // Member Operators

        //- Return element of UList
        inline scalar& operator[](const label);

        //- Return element of constant UList
        //  Note that the bool specialisation adds lazy evaluation so reading
        //  an out-of-range element returns false without any ill-effects
        inline const scalar& operator[](const label) const;

        //- Assignment of all entries to the given value
        void operator=(const scalar&);

        //- Assignment of all entries to zero
        void operator=(const zero);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CmptListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
