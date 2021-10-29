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
    Foam::energyImmersedBoundaryFvPatchField

Description


SourceFiles
    energyImmersedBoundaryFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef energyImmersedBoundaryFvPatchScalarField_H
#define energyImmersedBoundaryFvPatchScalarField_H

#include "immersedBoundaryFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class energyImmersedBoundaryFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

template<class BaseType>
class energyImmersedBoundaryFvPatchScalarField
:
    public BaseType
{
// Private Data

    //- Reference to global fvMesh
    const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("energy");


    // Constructors

        //- Construct from patch and internal field
        energyImmersedBoundaryFvPatchScalarField
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

#ifdef NoRepository
    #include "energyImmersedBoundaryFvPatchScalarField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
