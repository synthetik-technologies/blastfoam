/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023-2025
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
    Foam::errorTypes::field

Description
    Compute error in field estimates. Can use operations.

SourceFiles
    fieldError.C

\*---------------------------------------------------------------------------*/

#ifndef fieldError_H
#define fieldError_H

#include "errorType.H"
#include "runTimeSelectionTables.H"

#include "Function1Fwd.H"
#include "Function2Fwd.H"
#include "Function3Fwd.H"
#include "Function4Fwd.H"

#include "fvc.H"

namespace Foam
{
namespace errorTypes
{

/*---------------------------------------------------------------------------*\
                            Class field Declaration
\*---------------------------------------------------------------------------*/

class field
:
    public errorType
{

public:

    // Public enumerations
    enum Operation
    {
        NONE,
        LAPLACIAN,
        DDT,
        MAG_GRAD,
        FUNCTION
    };
    static const NamedEnum<Operation, 5> operationNames_;

    enum FieldReduction
    {
        MIN = Reduction::MIN,
        MIN_MAG_SQR = Reduction::MIN_MAG_SQR,
        MAX = Reduction::MAX,
        MAX_MAG_SQR = Reduction::MAX_MAG_SQR,
        AVERAGE = Reduction::AVERAGE,
        WEIGHTED_AVERAGE = Reduction::WEIGHTED_AVERAGE,

        // Replace "LATEST"
        PROBE = Reduction::LATEST
    };
    static const NamedEnum<FieldReduction, 7> fieldReductionTypeNames;

protected:
    // Protected data

        //- Name of the field
        const word fieldName_;

        //- Type of operation
        const Operation operation_;

        //- How the field is reduced
        const FieldReduction fieldReduction_;

        bool timeVarying_;
        labelList cmpts_;
        autoPtr<Function1<scalar>> func1_;
        autoPtr<Function2<scalar>> func2_;
        autoPtr<Function3<scalar>> func3_;
        autoPtr<Function4<scalar>> func4_;
        vector probeLocation_;


    // Protected member functions

        //- Perform operation
        tmp<volScalarField> performOperation(const volScalarField& f) const;

public:

    //- Runtime type information
    TypeName("field");

    // Constructor

        // Construct from dictionary
        field
        (
            const Time& runTime,
            const dictionary& dict,
            const word& region
        );

        // Construct from dictionary and a given operation
        field
        (
            const Time& runTime,
            const dictionary& dict,
            const word& region,
            const Operation op
        );

    //- Clone the object
    virtual autoPtr<errorType> clone() const
    {
        return autoPtr<errorType>(new field(*this));
    }


    //- Destructor
    virtual ~field();


    //Member functions

        //- Return the field name
        virtual const word& fieldName() const
        {
            return fieldName_;
        }

        //- Update the error
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace errorTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
