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
    Foam::equation

Description
    Scalar equation with a single input. Derived from both scalarEquation
    and multivariateEquation so that a scalar equation can be cast in either
    form.

\*---------------------------------------------------------------------------*/

#ifndef regEquation_H
#define regEquation_H

#include "regIOobject.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class regEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class BaseEquation>
class regEquation
:
    public regIOobject,
    public BaseEquation<Type>
{
public:
    // Public typedef
    typedef regEquation<Type, BaseEquation> RegEquationType;

protected:
    // Protected member data

        // Reference to the constructing object registry
        const objectRegistry& obr_;


public:

    //- Runtime type information
    TypeName("regEquation");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        RegEquationType,
        dictionary,
        (
            const objectRegistry& obr,
            const dictionary& dict
        ),
        (obr, dict)
    );

    // Constructors

        //- Construct given bounds
        regEquation(const objectRegistry& obr, const dictionary& dict);


    //- Destructor
    virtual ~regEquation();


    //- Selector
    static autoPtr<typename BaseEquation<Type>::BaseEquation> New
    (
        const word& type,
        const objectRegistry& obr,
        const dictionary& dict
    );

    virtual bool writeData(Ostream&) const
    {
        return true;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "regEquation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
