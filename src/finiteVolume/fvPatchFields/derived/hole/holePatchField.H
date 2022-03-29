/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::holePatchField

Description

Usage
    Base class for patchFields that have holes that restrict motion


SourceFiles
    holePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef holePatchField_H
#define holePatchField_H

#include "fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class holePatchField Declaration
\*---------------------------------------------------------------------------*/

class holePatchField
{
protected:
    // Protected Data

        //- Reference to the underlying patch
        const fvPatch& patch_;

        //- Name of the length field
        word LName_;

        //- Size of a hole (effective diameter)
        scalarField L_;

        //- Hole density
        scalarField holeDensity_;

        //- Fraction of surface that is open
        scalarField openFraction_;

        //- Mask (0 is wall, 1 is hole)
        scalarField mask_;


public:

    //- Runtime type information
    TypeName("hole");


    // Constructors

        //- Construct from patch
        holePatchField(const fvPatch&);

        //- Construct from patch and dictionary
        holePatchField(const fvPatch&, const dictionary&);

        //- Construct by mapping given holePatchField
        holePatchField
        (
            const holePatchField&,
            const fvPatch&,
            const fvPatchFieldMapper&
        );

    //- Destructor
    virtual ~holePatchField()
    {}


    // Member Functions

        // Access

            //- Return the mask
            const scalarField& mask() const
            {
                return mask_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            void rmap(const holePatchField&, const labelList&);

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            void updateCoeffs();

        //- Write
        void write(Ostream&) const;


    // Static Member functions

        // Remove this boundary from the given pressure gradient field
        static void removeBoundaryGradient
        (
            const volVectorField& U,
            volVectorField& gradP
        );

        // Remove this boundary from the given pressure gradient field
        static tmp<volVectorField> removeBoundaryGradient
        (
            const volVectorField& U,
            const tmp<volVectorField>& gradP
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
