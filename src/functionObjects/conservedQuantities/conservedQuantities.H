/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
28-10-2020 Synthetik Applied Technologies: | Calculate conservedQuantities
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
    Foam::functionObjects::conservedQuantities

Description
    Tracks the volume integrated quantities specified

    Example of function object specification:
    \verbatim
    conservedQuantities
    {
        type                conservedQuantities;
        fields              (rho rhoE rhoU);
    }
    \endverbatim

Usage
    \table
        Property          | Description               | Required | Default
        type              | type name: conservedQuantities | yes      |
        fields            | Tracked fields            | yes | no
    \endtable


See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    conservedQuantities.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_conservedQuantities_H
#define functionObjects_conservedQuantities_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{


/*---------------------------------------------------------------------------*\
                        Class conservedQuantities Declaration
\*---------------------------------------------------------------------------*/

class conservedQuantities
:
    public fvMeshFunctionObject
{
    // Private data

        //- Name of conserved field names
        wordList fields_;

        //- Scalar fields
        HashTable<dimensionedScalar> scalar0s_;
        HashTable<dimensionedScalar> scalars_;

        //- Vector fields
        HashTable<dimensionedVector> vector0s_;
        HashTable<dimensionedVector> vectors_;

        //- Symm tensor fields
        HashTable<dimensionedSymmTensor> symmTensor0s_;
        HashTable<dimensionedSymmTensor> symmTensors_;

        //- Spherical tensor fields
        HashTable<dimensionedSphericalTensor> sphericalTensor0s_;
        HashTable<dimensionedSphericalTensor> sphericalTensors_;

        //- Tensor fields
        HashTable<dimensionedTensor> tensor0s_;
        HashTable<dimensionedTensor> tensors_;


    // Private member functions 

        //- Functions to return the given table types
        //  these are later specified for all standard types
        template<class Type>
        const HashTable<dimensioned<Type>>& vals0() const
        {
            NotImplemented;
            return NullObjectRef<HashTable<dimensioned<Type>>>();
        }
        
        template<class Type>
        HashTable<dimensioned<Type>>& vals0()
        {
            NotImplemented;
            return const_cast<HashTable<dimensioned<Type>>&>
            (
                NullObjectRef<HashTable<dimensioned<Type>>>()
            );
        }

        template<class Type>
        const HashTable<dimensioned<Type>>& vals() const
        {
            NotImplemented;
            return NullObjectRef<HashTable<dimensioned<Type>>>();
        }

        template<class Type>
        HashTable<dimensioned<Type>>& vals()
        {
            NotImplemented;
            return NullObjectRef<HashTable<dimensioned<Type>>>();
        }


public:

    //- Runtime type information
    TypeName("conservedQuantities");


    // Constructors

        //- Construct from Time and dictionary
        conservedQuantities
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        conservedQuantities(const conservedQuantities&) = delete;


    //- Destructor
    virtual ~conservedQuantities();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Calculate the field averages
        virtual bool execute();

        //- Write the field averages
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const conservedQuantities&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
