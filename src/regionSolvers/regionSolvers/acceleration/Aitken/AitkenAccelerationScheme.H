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
    Foam::accelerationSchemes::Aitken

Description
    Aitken acceleration

SourceFiles
    AitkenAccelerationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef AitkenAccelerationScheme_H
#define AitkenAccelerationScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ResidualBaseAccelerationScheme.H"

namespace Foam
{
namespace accelerationSchemes
{

/*---------------------------------------------------------------------------*\
                           Class Aitken Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class Patch, class Mesh>
class Aitken
:
    public ResidualBase<Type, Patch, Mesh>
{
    //- Number of fixed relaxation steps
    label nFixed_;

    //- Initial relaxation factor
    scalar initRelaxFactor_;

    //- Maximum relaxation factor
    scalar maxRelaxFactor_;

    //- Current Aitken factor
    scalar aitkenFactor_;

public:

    TypeName("Aitken");

    // Constructor
    Aitken
    (
        GeometricField<Type, Patch, Mesh>& field,
        autoPtr<PatchFieldSelector<Type>> selector,
        const dictionary& dict
    );


    //- Destructor
    virtual ~Aitken();


    // Member functions

        //- Relax the field
        virtual void relax(const label iter);

        //- Read parameters
        virtual void read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace AccelerationSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AitkenAccelerationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
