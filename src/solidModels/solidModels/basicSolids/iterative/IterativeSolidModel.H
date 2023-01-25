/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    IterativeSolidModel

Description
    Virtual base class for solid mechanics models wit iterations

SourceFiles
    IterativeSolidModel.C
    IterativeSolidModelTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef IterativeSolidModel_H
#define IterativeSolidModel_H

#include "dynamicFvMesh.H"
#include "IOdictionary.H"
#include "scalarIOField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "thermalModel.H"
#include "mechanicalModel.H"
#include "nonLinearGeometry.H"
#include "uniformDimensionedFields.H"
#include "momentumStabilisation.H"
#include "setCellDisplacements.H"
#include "OFstream.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IterativeSolidModel Declaration
\*---------------------------------------------------------------------------*/

template<class SolidModel>
class IterativeSolidModel
:
    public SolidModel
{
    // Private data

        //- Solution standard tolerance
        scalar solutionTol_;

        //- Material law tolerance
        scalar materialTol_;

        //- Solution tighter tolerance
        scalar alternativeTol_;

        //- Write frequency for residuals information
        int infoFrequency_;

        //- Maximum number of momentum correctors
        int nCorr_;

        //- Minium number of momentum correctors
        int minCorr_;

        //- Number of times the maximum number of correctors was reached in the
        //  momentum equation
        int maxIterReached_;

        //- Residual file
        mutable autoPtr<OFstream> residualFilePtr_;


        //- Under-relaxation method
        //  It can be 'fixed', 'Aitken' or 'QuasiNewton'
        word relaxationMethod_;


        // Aitken fields

            //- Aitken's under-relaxation factor field
            autoPtr<volScalarField> aitkenAlpha_;

            //- Aitken's residual field
            autoPtr<volVectorField> aitkenResidual_;


        // Quasi-Newton fields

            //- Restart frequency for the Quasi-Newton method
            label QuasiNewtonRestartFreq_;

            //- Quasi-Newton method input vector field
            DynamicList<vectorField> QuasiNewtonV_;

            //- Quasi-Newton method output vector field
            DynamicList<vectorField> QuasiNewtonW_;

            //- Quasi-Newton method time indices for the input-ouput fields
            DynamicList<scalar> QuasiNewtonT_;

            //- Reference D field used for Quasi-Newton method
            autoPtr<volVectorField> DRef_;

            //- Reference D (without relaxation) field used for Quasi-Newton method
            autoPtr<volVectorField> unrelaxedDRef_;


protected:

        virtual void readDict();

        //- Apply fixed or Aitken's adaptive under-relaxation to the field
        void relaxField(volVectorField& D, const label iCorr);

        //- Check if the equation has converged
        template<class Type>
        bool converged
        (
            const label iCorr,
            const scalar solverPerfInitRes,
            const label solverPerfNIters,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const bool writeResiduals = true
        );

        //- Return the residual file output
        Ostream& residualFile() const
        {
            return residualFilePtr_();
        }

        //- Return solution standard tolerance
        scalar solutionTol() const
        {
            return solutionTol_;
        }

        //- Return solution tighter tolerance
        scalar alternativeTol() const
        {
            return alternativeTol_;
        }

        //- Return material law tolerance
        scalar materialTol() const
        {
            return materialTol_;
        }

        //- Return write frequency for residuals information
        int infoFrequency() const
        {
            return infoFrequency_;
        }

        //- Return maximum number of correctors
        int nCorr() const
        {
            return nCorr_;
        }

        //- Return maxIterReached
        int& maxIterReached()
        {
            return maxIterReached_;
        }


public:

    // Constructors

        //- Construct from components
        IterativeSolidModel
        (
            const word& type,
            dynamicFvMesh& mesh,
            const nonLinearGeometry::nonLinearType nonlinear,
            const bool incremental,
            const bool isSolid
        );


    // Destructor
    virtual ~IterativeSolidModel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
    #include "IterativeSolidModel.C"
    #include "IterativeSolidModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
