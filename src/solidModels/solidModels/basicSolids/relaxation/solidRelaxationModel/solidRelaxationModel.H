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
    solidRelaxationModel

Description
    Virtual base class for solid mechanics models wit iterations

SourceFiles
    solidRelaxationModel.C
    solidRelaxationModelTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef solidRelaxationModel_H
#define solidRelaxationModel_H

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
                           Class solidRelaxationModel Declaration
\*---------------------------------------------------------------------------*/

class solidRelaxationModel
{
protected:
    // Protected member data

        //- Reference to the solid model
        solidModel& solid_;


    // Protected member functions

        //- Read the dictionary
        virtual void readDict();

        //- Return the model dictionary
        const dictionary& coeffDict() const
        {
            return solid_.solidModelDict().optionalSubDict(type() + "Coeffs");
        }


public:

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidRelaxationModel,
            dictionary,
            (
                solidModel& solid
            ),
            (solid)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidRelaxationModel,
            dictionary,
            (
                solidModel& solid
            ),
            (solid)
        );

    // Constructors

        //- Construct from components
        solidRelaxationModel(solidModel& solid);


    // Selectors

        //- Select constructed from time
        static autoPtr<solidRelaxationModel> New(solidModel& solid);


    // Destructor
    virtual ~solidRelaxationModel();


    //- Apply under-relaxation to the field
    virtual void relaxField(volVectorField& D, const label iCorr) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
