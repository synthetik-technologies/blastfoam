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
    Foam::levelSetModel

Description
    Class to handle creation of level set and redistancing

SourceFiles
    levelSetModel.C

\*---------------------------------------------------------------------------*/

#ifndef levelSetModel_H
#define levelSetModel_H

#include "volFields.H"
#include "isoSurface.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class levelSetModel Declaration
\*---------------------------------------------------------------------------*/

class levelSetModel
{
public:

    enum levelSetFunc
    {
        TANH,
        EXP
    };

    //- Named enumeration of level set functions
    static const NamedEnum<levelSetFunc, 2> levelSetFuncNames_;

private:
    // Private Data

        //- Constant reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the volume fraction field
        const volScalarField& alpha_;

        //- Level set function
        volScalarField levelSet_;

        //- Thickness of interface
        dimensionedScalar epsilon_;

        //- Power used to smooth the levelSet function for more diffuse
        //  interface. This is used to compute surface normals
        scalar smoothPow_;

        //- Filtering method for contouring
        isoSurface::filterType filterType_;

        //- Method for computing level set function
        levelSetFunc lsFunc_;

        //- Dictionary used to construct distributedTriSurfaceMesh
        dictionary triMeshDict_;



    // Private member functions

        //- Compute the level set function
        tmp<volScalarField> calcLevelSet(const volScalarField&) const;

public:

    //- Runtime type information
    TypeName("levelSet");


    // Constructors

        //- Construct from volume fraction and dictionary
        levelSetModel
        (
            const volScalarField& alpha,
            const dictionary& dict
        );

    //- Destructor
    virtual ~levelSetModel();


    // Member Functions

        //- Recompute the levelSet function
        void redistance();

        //- Interface indictator function
        tmp<volScalarField> H() const;

        //- Return the smoothed levelSet function
        tmp<volScalarField> psi() const;

        //- Return the gradient of the smoothed level set function
        tmp<volVectorField> gradLevelSet() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
