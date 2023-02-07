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
    Foam::fieldAccelerationScheme

Description
    Base class for acceleration

SourceFiles
    fieldAccelerationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef fieldAccelerationScheme_H
#define fieldAccelerationScheme_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "accelerationScheme.H"
#include "regionSolverList.H"
#include "PtrList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fieldAccelerationScheme Declaration
\*---------------------------------------------------------------------------*/

class fieldAccelerationScheme
:
    public PtrList<accelerationScheme>
{
protected:
// Protected data

    // Field name
    const word& fieldName_;

    //- Reference to the mesh
    const fvMesh& mesh_;

    //- Tolerances
    scalar tolerance_;
    scalar relTol_;

    //- Read controls
    bool readControls(const dictionary& dict);

public:

    TypeName("fieldAccelerationScheme");

    // Constructor
    fieldAccelerationScheme
    (
        const word& field,
        const fvMesh& mesh,
        const labelList& patches,
        const dictionary& dict
    );

    //- Destructor
    virtual ~fieldAccelerationScheme();


    // Member functions

        //- Update residuals
        void updateError();

        //- Relax the field
        void relax(const label iter);

        //- Clear residuals
        void clear();

        // Return the absolute error
        scalar error() const;

        //- Return the relative error
        scalar relError() const;

        Convergence convergence() const
        {
            if (tolerance_ < 0)
            {
                return UNKNOWN_CONVERGENCE;
            }
            else if (error() < small)
            {
                return FULL_CONVERGENCE;
            }
            else if (error() < tolerance_ || relError() < relTol_)
            {
                return CONVERGED;
            }
            return NOT_CONVERGED;
        }

        //- Read parameters
        void read(const dictionary& dict);

        //- Print current errors
        void print(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
