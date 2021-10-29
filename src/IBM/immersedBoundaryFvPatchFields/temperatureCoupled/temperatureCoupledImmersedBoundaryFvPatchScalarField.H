/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
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
    Foam::temperatureCoupledImmersedBoundaryFvPatchField

Description


SourceFiles
    temperatureCoupledImmersedBoundaryFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureCoupledImmersedBoundaryFvPatchScalarField_H
#define temperatureCoupledImmersedBoundaryFvPatchScalarField_H

#include "immersedBoundaryFvPatchFields.H"
#include "immersedMeshMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class temperatureCoupledImmersedBoundaryFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class temperatureCoupledImmersedBoundaryFvPatchScalarField
:
    public immersedBoundaryScalarPatchField
{
// Private Data

    //- Reference to global fvMesh
    const fvMesh& mesh_;

    //- Reference to immersedMesh
    const fvMesh& immersedMesh_;

    //- Reference to immersed temperature field
    const volScalarField& immersedT_;

    //- Mapped used to map to and from the immersed mesh
    const immersedMeshMapper& mapper_;

    //- Phase name
    word phaseName_;

    //- Return kappa field
    tmp<scalarField> kappaImmersed() const;
    tmp<scalarField> kappa() const;


public:

    //- Runtime type information
    TypeName("temperatureCoupled");


    // Constructors

        //- Construct from patch and internal field
        temperatureCoupledImmersedBoundaryFvPatchScalarField
        (
            volScalarField&,
            const dictionary&,
            const immersedBoundaryObject&
        );


    // Member Functions


        // Evaluation functions

            //- Update size and values
            virtual void updateCoeffs() const;

            //- Add forcing to the given field
            virtual void setValues();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
