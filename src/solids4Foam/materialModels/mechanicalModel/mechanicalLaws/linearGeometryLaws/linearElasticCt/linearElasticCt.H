/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    linearElasticCt

Description
    Hookean linear elastic mechanical law.
    Sets Youngs Modulus of each cell from CT Data based on an empircal law.

SourceFiles
    linearElasticCt.C

Author
    Karen Fitzgerald, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef linearElasticCt_H
#define linearElasticCt_H

#include "mechanicalLaw.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElasticCt Declaration
\*---------------------------------------------------------------------------*/

class linearElasticCt
:
    public mechanicalLaw
{
    // Private data

        //- Young's modulus field
        volScalarField E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- Shear modulus field
        volScalarField mu_;

        //- 1st Lame parameter field
        volScalarField lambda_;

        //- Shear modulus surface field
        surfaceScalarField muf_;

        //- 1st Lame parameter surface field
        surfaceScalarField lambdaf_;

        //- Rotation mesh coordinates before looking CT values
        const Switch useRotationMatrix_;

        //- Rotation matrix to rotate current mesh back to CT reference frame
        tensor rotationMatrix_;

        //- Centre of rotation for transformation
        vector centreOfRotation_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        linearElasticCt(const linearElasticCt&);

        //- Disallow default bitwise assignment
        void operator=(const linearElasticCt&);

        //- Set E field from CT Data
        void setYoungsModulusFromCt();


public:

    //- Runtime type information
    TypeName("linearElasticCt");

    // Static data members


    // Constructors

        //- Construct from dictionary
        linearElasticCt
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~linearElasticCt();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> impK() const;

        //- Calculate stress volField
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate stress surfaceField
        virtual void correct(surfaceSymmTensorField& sigma);


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
