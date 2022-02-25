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
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class levelSetModel Declaration
\*---------------------------------------------------------------------------*/

class levelSetModel
{
public:

    enum class levelSetFunc
    {
        TANH,
        EXP
    };

    enum class truncation
    {
        NONE,
        CUTOFF,
        TANH
    };

    //- Named enumeration of level set functions
    static const NamedEnum<levelSetFunc, 2> levelSetFuncNames_;

    //- Named enumeration of level set truncation methods
    static const NamedEnum<truncation, 3> truncationNames_;

private:
    // Private Data

        //- Constant reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the volume fraction field
        volScalarField& alpha_;

        //- Level set function
        volScalarField levelSet_;

        //- Optional truncated levelSet
        autoPtr<volScalarField> tlevelSet_;

        //- Interface normal direction
        surfaceScalarField nHatf_;

        //- Interface curvature
        volScalarField K_;

        //- Thickness of interface
        dimensionedScalar epsilon0_;

        //- True interface thickness (dx*epsilon0)
        volScalarField epsilon_;

        //- Power used to smooth the levelSet function for more diffuse
        //  interface. This is used to compute surface normals
        scalar smoothPow_;

        //- Use distributed triSurface
        bool useDistributed_;

        //- Filtering method for contouring
        isoSurface::filterType filterType_;

        //- Method for computing level set function
        levelSetFunc lsFunc_;

        //- Method for truncating the level set function
        truncation truncation_;

        //- Truncation cut off value
        scalar cutOff_;

        //- Dictionary used to construct distributedTriSurfaceMesh
        mutable dictionary triMeshDict_;

        //- Stored isoSurface
        mutable autoPtr<isoSurface> contourPtr_;


public:

    //- Runtime type information
    TypeName("levelSet");


    // Constructors

        //- Construct from volume fraction and dictionary
        levelSetModel
        (
            volScalarField& alpha,
            const dictionary& dict
        );

    //- Destructor
    virtual ~levelSetModel();


    // Member Functions

        // Access

            //- Return the levelSet function
            const volScalarField& levelSet() const
            {
                return levelSet_;
            }

            //- Access the levelSet function
            volScalarField& levelSet()
            {
                return levelSet_;
            }

            //- Return the truncated levelSet function
            const volScalarField& tlevelSet() const
            {
                return tlevelSet_.valid() ? tlevelSet_() : levelSet_;
            }

            //- Return the interface thickness parameter
            const volScalarField& epsilon() const
            {
                return epsilon_;
            }

            //- Return the interface normal direction
            const surfaceScalarField& nHatf() const
            {
                return nHatf_;
            }

            //- Return the interface curvature
            const volScalarField& K() const
            {
                return K_;
            }

            //- Return the isoSurface
            const isoSurface& contour() const;


        //- Update and optionally calculate distance
        void correct(const bool redistance = false);

        //- Compute the level set function
        tmp<volScalarField> calcLevelSet(const volScalarField&) const;

        //- Compute the level set function
        tmp<volScalarField> calcDistance(const volScalarField&) const;

        //- Return the calculated volume fraction (0-1)
        tmp<volScalarField> alpha() const;

        //- Calculate and return the distance from the interface
        tmp<volScalarField> distance() const;

        //- Interface indictator function
        tmp<volScalarField> nearInterface() const;

        //- Interface cells
        labelList interfaceCells() const;

        //- Return the smoothed levelSet function
        tmp<volScalarField> psi() const;

        //- Return the gradient of the smoothed level set function
        tmp<volVectorField> gradLevelSet() const;

        //- Return the gradient of the volume fraction
        tmp<volVectorField> gradAlpha() const;

        //- Return interface normal
        tmp<volVectorField> nHat() const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
