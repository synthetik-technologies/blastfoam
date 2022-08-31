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
#include "searchableSurface.H"

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

        //- Heavy-side function
        volScalarField H_;

        //- Interface normal direction
        surfaceScalarField nHatf_;

        //- Interface curvature
        volScalarField K_;

        //- Thickness of interface
        dimensionedScalar epsilon0_;

        //- True interface thickness (dx*epsilon0)
        volScalarField epsilon_;

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


public:

    //- Runtime type information
    TypeName("levelSet");


    // Constructors

        //- Construct from volume fraction and dictionary
        levelSetModel
        (
            volScalarField& alpha,
            const dictionary& dict,
            const bool mustRead = true
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

            //- Return the heavy-side function
            const volScalarField& H() const
            {
                return H_;
            }

            //- Access the heavy-side function
            volScalarField& H()
            {
                return H_;
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


        //- Update and optionally calculate distance
        void correct(const bool updateH = true);

        //- Update the epsilon field
        void updateEpsilon();

        //- Compute the level set function
        tmp<volScalarField> calcH(const volScalarField& ls) const;

        //- Compute the distance function using the inverse levelSet function
        tmp<volScalarField> calcLevelSet(const volScalarField& H) const;

        // Calculate the levelSet from a list of searchableSurfaces
        tmp<volScalarField> calcLevelSet
        (
            const volScalarField& alpha,
            const UPtrList<searchableSurface>& regions
        ) const;

        //- Calculate and return the distance from the interface
        tmp<volScalarField> calcLevelSet
        (
            const volScalarField& isoField,
            const scalar isoValue
        ) const;

        //- Recompute distance from levelSet 0
        void redistance();

        //- Return the calculated volume fraction (0-1)
        tmp<volScalarField> alpha() const;

        //- Interface indictator function
        tmp<volScalarField> nearInterface() const;

        //- Return interface normal
        tmp<volVectorField> nHat() const;

        //- Return the gradient of the volume fraction
        tmp<volVectorField> gradAlpha() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
