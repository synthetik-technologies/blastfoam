/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2020
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
    Foam::errorEstimator

Description
    Base class used to estimate error within a cell/across faces

SourceFiles
    errorEstimator.C
    newErrorEstimator.C

\*---------------------------------------------------------------------------*/

#ifndef errorEstimator_H
#define errorEstimator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class errorEstimator Declaration
\*---------------------------------------------------------------------------*/

class errorEstimator
{
protected:
    // Protected data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Name of error estimator
        word name_;

        //- Reference to error field
        volScalarField& error_;

        //- Lower refine level
        scalar lowerRefine_;

        //- Lower unrefine level
        scalar lowerUnrefine_;

        //- Upper refine level
        scalar upperRefine_;

        //- Upper unrefine level
        scalar upperUnrefine_;

        //- Create or lookup error field
        //  Useful when multiple error estimators are used
        volScalarField& lookupOrConstruct(const fvMesh& mesh) const;

        //- Normalize error (-1 = unrefine, 0 = do nothing, 1 = refine)
        void normalize(volScalarField& error);

        //- Lookup and return the scalar value of the field
        void getFieldValue(const word& name, volScalarField& f) const;

        //- Lookup and return the scalar value of the field
        template<class Type>
        void getFieldValueType(const word& name, volScalarField& f) const;

public:

    //- Runtime type information
    TypeName("errorEstimator");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            errorEstimator,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const word& name
            ),
            (mesh, dict, name)
        );

    // Constructor
    errorEstimator
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& name
    );


    //- Destructor
    virtual ~errorEstimator();


    // Selectors

        static autoPtr<errorEstimator> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& name = word::null
        );


    // Member Functions

        //- Update
        virtual void update() = 0;

        //- Return constant reference to error field
        const volScalarField& error() const
        {
            return error_;
        }

        //- Return non constant reference to error field
        volScalarField& error()
        {
            return error_;
        }

        //- Read level values
        virtual void read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "errorEstimatorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
