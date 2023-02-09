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
    Foam::accelerationSchemes::NoAcceleration

Description
    No acceleration

SourceFiles
    NoAccelerationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef NoAccelerationScheme_H
#define NoAccelerationScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AccelerationSchemeBase.H"

namespace Foam
{
namespace accelerationSchemes
{

/*---------------------------------------------------------------------------*\
                           Class NoAcceleration Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class Patch, class Mesh>
class NoAcceleration
:
    public AccelerationSchemeBase<Type, Patch, Mesh>
{
public:

    TypeName("none");

    // Constructor
    NoAcceleration
    (
        GeometricField<Type, Patch, Mesh>& field,
        autoPtr<PatchFieldSelector<Type>> selector,
        const dictionary& dict
    );


    //- Destructor
    virtual ~NoAcceleration();


    // Member functions

        //- Relax the field
        virtual void relax(const label iter);

        //- Read parameters
        virtual void read(const dictionary& dict)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace accelerationSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoAccelerationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
