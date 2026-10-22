/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
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
    Foam::errorEstimators::linearCombination

Description
    Error is a linear combination of scalar fields, optionally divided by a
    second linear combination:

        error = sum_i a_i f_i / sum_j b_j f_j

    Cells with lowerRefineLevel < error < upperRefineLevel are refined, cells
    with error < unrefineLevel or error > upperUnrefineLevel are coarsened.
    With the element mass fractions of the species, e.g. the ratio of the
    carbon to the carbon and nitrogen of the gas, it marks a mixing layer
    between detonation products and air independently of the chemistry.
    Fields that are not (yet) registered are skipped.

Usage
    \verbatim
    mixingLayer
    {
        errorEstimator      linearCombination;
        fields              { CO.gas 0.4288; CO2.gas 0.2729; }
        denominator         { CO.gas 0.4288; CO2.gas 0.2729; N2.gas 1; }
        small               1e-10;      // Floor of the denominator
        lowerRefineLevel    0.01;
        upperRefineLevel    0.44;
        unrefineLevel       0.005;
        upperUnrefineLevel  0.445;
    }
    \endverbatim

SourceFiles
    linearCombinationErrorEstimator.C

\*---------------------------------------------------------------------------*/

#ifndef linearCombinationErrorEstimator_H
#define linearCombinationErrorEstimator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "errorEstimator.H"

namespace Foam
{
namespace errorEstimators
{

/*---------------------------------------------------------------------------*\
                      Class linearCombination Declaration
\*---------------------------------------------------------------------------*/

class linearCombination
:
    public errorEstimator
{
// Private data

        //- Names and coefficients of the numerator fields
        wordList fields_;
        scalarList coeffs_;

        //- Names and coefficients of the denominator fields (optional)
        wordList denFields_;
        scalarList denCoeffs_;

        //- Floor of the denominator
        scalar small_;


    // Private Member Functions

        //- Read the field names and coefficients of a dictionary
        static void readCoeffs
        (
            const dictionary& dict,
            wordList& fields,
            scalarList& coeffs
        );

        //- Sum of the registered fields with their coefficients
        tmp<volScalarField> sum
        (
            const wordList& fields,
            const scalarList& coeffs
        ) const;


public:

    //- Runtime type information
    TypeName("linearCombination");

    // Constructor
    linearCombination
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& name
    );


    //- Destructor
    virtual ~linearCombination();


    // Member Functions

        //- Update error
        virtual void update(const bool scale = true);

        //- Read parameters
        virtual void read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace errorEstimators
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
