/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::plasticSolver

Description
    Class for solving plastic deformation

SourceFiles
    plasticSolver.C

\*---------------------------------------------------------------------------*/

#ifndef plasticSolver_H
#define plasticSolver_H

#include "materialModel.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class plasticSolver Declaration
\*---------------------------------------------------------------------------*/

class plasticSolver
{
protected:
// Protected data

    static const scalar sqrt2By3;

    //- Maximum allowed error in the plastic strain integration
    const scalar maxDeltaErr_;

    //- Tolerance for Newton loop
    const scalar tolerance_;

    //- Maximum number of iterations for Newton loop
    const label maxIter_;

    //- finiteDiff is the delta for finite difference differentiation
    const scalar finiteDifference_;


public:

    //- Runtime type information
    ClassName("plasticSolver");

    // Constructors

        //- Construct from order, size and element type
        plasticSolver(const dictionary& dict);

    //- Destructor
    virtual ~plasticSolver();

    struct yieldStressModel
    :
        public tmp<yieldStressModel>::refCount
    {
        yieldStressModel()
        {}

        virtual ~yieldStressModel()
        {}

        //- Return the current yield stress
        virtual void update(const label, const label)
        {}

        //- Return the current yield stress
        virtual scalar currentSigmaY
        (
            const scalar epsilonPEq0,
            const scalar curentEpsilonPEq,
            const scalar J
        ) const = 0;
    };


    // Member functions

        // Evaluation

            //- Return a yield stress model
            virtual tmp<yieldStressModel> yieldStress
            (
                const label elemi,
                const label rulei
            ) const = 0;

            //- Evaulate current value of the yield function
            virtual scalar fY
            (
                const yieldStressModel& ys,
                const scalar epsilonPEq0,
                const scalar magSTrial,
                const scalar DLambda,
                const scalar muBar,
                const scalar J
            ) const;

            //- Iteratively calculate plastic multiplier increment (DLambda)
            //  and current yield stress using Newton's method
            virtual void newtonLoop
            (
                const yieldStressModel& ys,
                scalar& DLambda,
                scalar& curentSigmaY,
                const scalar epsilonPEq0,
                const scalar magSTrial,
                const scalar muBar,
                const scalar maxMagDEpsilon,
                const scalar J
            ) const;

            //- Update plasticN, DLambda, DSigmaY and sigmaY
            virtual void updatePlasticity
            (
                const yieldStressModel& ys,
                symmTensor& plasticN,
                scalar& DLambda,
                scalar& sigmaY,
                const scalar sigmaY0,
                const scalar fTrial,
                const symmTensor& sTrial,
                const scalar epsilonPEq0,
                const scalar muBar,
                const scalar maxMagDEpsilon,
                const scalar J
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
