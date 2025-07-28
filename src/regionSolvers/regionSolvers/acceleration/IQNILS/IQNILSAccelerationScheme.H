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
    Foam::accelerationSchemes::IQNILS

Description
    IQNILS acceleration

SourceFiles
    IQNILSAccelerationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef IQNILSAccelerationScheme_H
#define IQNILSAccelerationScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "QNBaseAccelerationScheme.H"
#include "NamedEnum.H"

namespace Foam
{
namespace accelerationSchemes
{

/*---------------------------------------------------------------------------*\
                           Class IQNILS Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class Patch, class Mesh>
class IQNILS
:
    public QNBase<Type, Patch, Mesh>
{
private:
// Private Data

    //- Initial relaxation factor
    scalar initRelaxFactor_;

    //- Number of fixed relaxation steps
    label nFixed_;

    //- Reference field (unrelaxed)
    Field<Type> refField_;

    //- QR decomposition method
    word qr_;


    // Private member function

        //- Return V reference residuals
        virtual const Field<Type>& VRef() const
        {
            return this->prevResiduals_;
        }

        //- Return W reference field
        virtual const Field<Type>& WRef() const
        {
            return refField_;
        }

        //- QR decomposition using modified Gram-Schmidt
        void GramSchmidt_QRdecomp(Field<Type>& pFieldNew);


public:

    TypeName("IQNILS");

    // Constructor
    IQNILS
    (
        GeometricField<Type, Patch, Mesh>& field,
        autoPtr<PatchFieldSelector<Type>> selector,
        const dictionary& dict
    );


    //- Destructor
    virtual ~IQNILS();


    // Member functions

        //- Relax the field
        virtual void relax(const label iter);

        //- Read parameters
        virtual void read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace accelerationSchemes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IQNILSAccelerationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
