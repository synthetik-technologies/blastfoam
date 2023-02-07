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
    Foam::accelerationSchemes::QNBase

Description
    QNBase acceleration

SourceFiles
    QNBaseAccelerationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef QNBaseAccelerationScheme_H
#define QNBaseAccelerationScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ResidualBaseAccelerationScheme.H"

namespace Foam
{
namespace accelerationSchemes
{

/*---------------------------------------------------------------------------*\
                           Class QNBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class Patch, class Mesh>
class QNBase
:
    public ResidualBase<Type, Patch, Mesh>
{
protected:

    label maxCouplingIter_;
    label nCouplingTimes_;

    DynamicList<Field<Type>> Vs_;
    DynamicList<Field<Type>> Ws_;
    DynamicList<label> times_;

    //- Update matrices
    void updateVW(const label iter);

    //- Return V reference residuals
    virtual const Field<Type>& VRef() const = 0;

    //- Return W reference field
    virtual const Field<Type>& WRef() const = 0;

    //- Filter
    void filter();

public:

    // Constructor
    QNBase
    (
        const word& type,
        GeometricField<Type, Patch, Mesh>& field,
        const label patchi,
        const dictionary& dict
    );


    //- Destructor
    virtual ~QNBase();


    // Member functions

        //- Clear residuals
        virtual void clear();

        //- Read parameters
        virtual void read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace accelerationScehemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "QNBaseAccelerationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
