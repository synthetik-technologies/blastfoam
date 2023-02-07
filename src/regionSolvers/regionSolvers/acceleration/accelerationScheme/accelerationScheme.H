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
    Foam::accelerationScheme

Description
    Base class for acceleration

SourceFiles
    accelerationScheme.C
    accelerationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef accelerationScheme_H
#define accelerationScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "valuePointPatchField.H"
#include "runTimeSelectionTables.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class accelerationScheme Declaration
\*---------------------------------------------------------------------------*/

class accelerationScheme
{
protected:
// Protected data

        const word& type_;

        //- Patches to relax
        const label patchi_;

        //- Initial error
        scalar initialError_;

        //- Current residual
        scalar error_;

        //- Sum magSqr
        template<class Type>
        scalar sumDotDot(const UList<Type>& f) const
        {
            scalar sumf = Zero;
            forAll(f, i)
            {
                sumf += magSqr(f[i]);
            }
            return returnReduce(sumf, sumOp<scalar>());
        }

        template<class Type>
        scalar sumDotDot(const UList<Type>& f1, const UList<Type>& f2) const
        {
            scalar sumf12 = Zero;
            forAll(f1, i)
            {
                sumf12 += f1[i] && f2[i];
            }
            return returnReduce(sumf12, sumOp<scalar>());
        }

        template<class Type>
        scalar dotDot(const Type& f1) const
        {
            return magSqr(f1);
        }

        template<class Type>
        scalar dotDot(const Type& f1, const Type& f2) const
        {
            return f1 && f2;
        }


public:

    TypeName("acceleration");

    // Constructor
    accelerationScheme
    (
        const word& type,
        const label patchi,
        const dictionary& dict
    );


    // Selector
        template<class Type, template<class> class Patch, class Mesh>
        static autoPtr<accelerationScheme> New
        (
            GeometricField<Type, Patch, Mesh>& field,
            const label patchi,
            const dictionary& dict
        );


    //- Destructor
    virtual ~accelerationScheme();


    // Member functions

        //- Update residuals
        virtual void updateError() = 0;

        //- Relax the field
        virtual void relax(const label iter) = 0;

        //- Clear residuals
        virtual void clear();

        // Return the absolute error
        inline scalar error() const
        {
            return mag(error_);
        }

        // Return the initial absolute error
        inline scalar initError() const
        {
            return mag(initialError_);
        }

        //- Return the relative error
        inline scalar relError() const
        {
            return error_/(initialError_ + small);
        }

        //- Read parameters
        virtual void read(const dictionary& dict) = 0;
};

template<>
Foam::scalar Foam::accelerationScheme::sumDotDot
(
    const UList<scalar>& f
) const;
template<>
Foam::scalar Foam::accelerationScheme::sumDotDot
(
    const UList<scalar>& f1,
    const UList<scalar>& f2
) const;
template<>
Foam::scalar Foam::accelerationScheme::dotDot
(
    const scalar& f
) const;
template<>
Foam::scalar Foam::accelerationScheme::dotDot
(
    const scalar& f1,
    const scalar& f2
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
