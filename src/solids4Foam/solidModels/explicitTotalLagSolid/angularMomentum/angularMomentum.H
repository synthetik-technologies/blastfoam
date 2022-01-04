/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::angularMomentum

Description
    Angular momentum preservation algorithm

SourceFiles
    angularMomentum.C

\*---------------------------------------------------------------------------*/

#ifndef angularMomentum_H
#define angularMomentum_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class angularMomentum
{
    // Private data

        //- Mesh
        const fvMesh& mesh_;

        //- Material density
        const volScalarField& rho_;

        //- Radial axes
        vector radialAxes_;

        //- Rotation axes
        vector rotationAxes_;

        //- Valid directions (geometricD)
        vector validD_;

        //- Ignored directions (inverse of geometricD)
        vector invalidD_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        angularMomentum(const angularMomentum&);

        //- Disallow default bitwise assignment
        void operator=(const angularMomentum&);

public:

    //- Runtime type information
    TypeName("angularMomentum");


    // Constructors

        //- Construct from mesh and dictionary
        explicit angularMomentum
        (
            const fvMesh&,
            const dictionary&
        );


    //- Destructor
    virtual ~angularMomentum();


    // Member functions

        // Edit

            //- Angular momentum conservation
            void AMconservation
            (
                GeometricField<vector, fvPatchField, volMesh>&,
                GeometricField<vector, fvPatchField, volMesh>&,
                const GeometricField<vector, fvPatchField, volMesh>&,
                const label
            ) const;

        // Print

            //- Print global momentum
            void printGlobalMomentum
            (
                const GeometricField<vector, fvPatchField, volMesh>&,
                const GeometricField<vector, fvPatchField, volMesh>&
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
