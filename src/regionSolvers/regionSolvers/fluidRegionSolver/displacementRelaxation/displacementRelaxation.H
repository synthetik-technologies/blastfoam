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
    Foam::displacementRelaxation

Description
    Base class for relaxation

SourceFiles
    displacementRelaxation.C
    displacementRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef displacementRelaxation_H
#define displacementRelaxation_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dynamicFvMesh.H"
#include "pointFields.H"
#include "globalPolyBoundaryMesh.H"
#include "runTimeSelectionTables.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class displacementRelaxation Declaration
\*---------------------------------------------------------------------------*/

class displacementRelaxation
{
protected:
// Protected data

        //- fvMesh
        const fvMesh& mesh_;

        //- Indices of coupled patches
        labelList coupledPatches_;

        //- Initial error
        scalar initialError_;

        //- Current residual
        scalar error_;

        //- Return the coefficient dictionary
        const dictionary& coeffDict(const dictionary& dict) const
        {
            return dict.optionalSubDict(this->type() + "Coeffs");
        }


public:

    TypeName("displacementRelaxation");

     // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            displacementRelaxation,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );

    // Constructor
    displacementRelaxation(const fvMesh& mesh, const dictionary& dict);


    // Selector
        static autoPtr<displacementRelaxation> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~displacementRelaxation();


    // Member functions

        //- Update residuals
        virtual void updateError(const pointVectorField&);

        //- Relax the field
        virtual void relax(const label iter, pointVectorField&) = 0;

        //- Clear residuals
        virtual void clear();

        // Return the absolute error
        inline scalar error() const
        {
            return error_;
        }

        //- Return the relative error
        inline scalar relError() const
        {
            return error_/(initialError_ + small);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
