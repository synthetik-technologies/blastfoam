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
    Foam::accelerationSchemes::ResidualBase

Description
    Residual based acceleration

SourceFiles
    ResidualBaseAccelerationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef ResidualBaseAccelerationScheme_H
#define ResidualBaseAccelerationScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AccelerationSchemeBase.H"

namespace Foam
{
namespace accelerationSchemes
{

/*---------------------------------------------------------------------------*\
                           Class ResidualBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class Patch, class Mesh>
class ResidualBase
:
    public AccelerationSchemeBase<Type, Patch, Mesh>
{
protected:

    Field<Type> residuals_;
    Field<Type> prevResiduals_;


public:

    // Constructor
    ResidualBase
    (
        const word& type,
        GeometricField<Type, Patch, Mesh>& field,
        autoPtr<PatchFieldSelector<Type>> selector,
        const dictionary& dict
    );


    //- Destructor
    virtual ~ResidualBase();


    // Member functions

        //- Update residuals
        virtual void updateError();

        //- Clear residuals
        virtual void clear(const bool);

        //- Read parameters
        virtual void read(const dictionary& dict)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace accelerationScehemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ResidualBaseAccelerationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
