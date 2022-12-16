/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    unsSolidModel

Description
    Base class for unstructured solid models

SourceFiles
    unsSolidModel.C

\*---------------------------------------------------------------------------*/

#ifndef unsSolidModel_H
#define unsSolidModel_H

#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class unsSolidModel Declaration
\*---------------------------------------------------------------------------*/

class unsSolidModel
:
    public solidModel
{
    // Protected data

        //- Cauchy (true) stress surface field
        surfaceSymmTensorField sigmaf_;

        //- Surface displacement gradient
        surfaceTensorField gradDf_;

        //- Surface displacement increment gradient
        surfaceTensorField gradDDf_;


public:

    // Constructors

        //- Construct from components
        unsSolidModel
        (
            const word& type,
            dynamicFvMesh& mesh,
            const nonLinearGeometry::nonLinearType nonlinear,
            const bool incremental,
            const bool isSolid = true
        );

    // Destructor

        virtual ~unsSolidModel()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether gradDf or gradDDf
            //  is the primary solution variable
            virtual const surfaceTensorField& solutionGradDf() const = 0;

            //- Return face gradient of displacement
            inline const surfaceTensorField& gradDf() const
            {
                return gradDf_;
            }

            //- Access face gradient of displacement
            inline surfaceTensorField& gradDf()
            {
                return gradDf_;
            }

            //- Return face gradient of delta displacement
            inline const surfaceTensorField& gradDDf() const
            {
                return gradDDf_;
            }

            //- Access face gradient of delta displacement
            inline surfaceTensorField& gradDDf()
            {
                return gradDDf_;
            }

            //- Return face stress
            inline const surfaceSymmTensorField& sigmaf() const
            {
                return sigmaf_;
            }

            //- Access face stress
            inline surfaceSymmTensorField& sigmaf()
            {
                return sigmaf_;
            }

            //- Access face stress
            using solidModel::sigma;
            inline const symmTensorField& sigma(const fvPatch& patch) const
            {
                return sigmaf_.boundaryField()[patch.index()];
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
