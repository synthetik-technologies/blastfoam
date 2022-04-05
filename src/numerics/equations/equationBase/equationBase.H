/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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
    Foam::equationBase

Description
    Base class for all equation types

\*---------------------------------------------------------------------------*/

#ifndef equationBase_H
#define equationBase_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class equationBase Declaration
\*---------------------------------------------------------------------------*/

class equationBase
{
protected:

    // Protected data

        //- Number in input variables
        word name_;


public:

    // Constructors

        //- Construct null
        equationBase();

        //- Construct null
        equationBase(const string& name);


    //- Destructor
    virtual ~equationBase();


    // Member functions

        //- Return the name
        virtual const word& name() const
        {
            return name_;
        }

        //- Access the name
        virtual word& name()
        {
            return name_;
        }

        //- Return the number of derivatives
        virtual label nDerivatives() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
