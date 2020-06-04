/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::kineticTheoryModels::frictionalStressModels::Princeton

Description

    Reference:
    \verbatim
        Agrawal, K., Loezos, P. N., Syamlal, M. & Sundaresan, S. (2001).
        "The role of meso-scale structures in rapid gas–solid flows",
        Journal of Fluid Mechanics, 445, 151-185.
    \endverbatim

SourceFiles
    PrincetonFrictionalStress.C

\*---------------------------------------------------------------------------*/

#ifndef PrincetonFrictionalStress_H
#define PrincetonFrictionalStress_H

#include "frictionalStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{

/*---------------------------------------------------------------------------*\
               Class Princeton Declaration
\*---------------------------------------------------------------------------*/

class Princeton
:
    public frictionalStressModel
{
    // Private data

        dictionary coeffDict_;

        //- Material constant for frictional normal stress
        dimensionedScalar Fr_;

        //- Material constant for frictional normal stress
        dimensionedScalar eta_;

        //- Material constant for frictional normal stress
        dimensionedScalar p_;

        //- Angle of internal friction
        dimensionedScalar phi_;

        //- Lower limit for (alphaMax - alpha1)
        dimensionedScalar alphaDeltaMin_;

        //- Lower limit when friction is "turned on"
        dimensionedScalar alphaMinFriction_;

        //- Lower limit when Schaeffer friction is "turned on"
        dimensionedScalar alphaMinFrictionByAlphap_;


public:

    //- Runtime type information
    TypeName("Princeton");


    // Constructors

        //- Construct from components
        Princeton(const dictionary& dict);


    //- Destructor
    virtual ~Princeton();


    // Member functions

        virtual tmp<volScalarField> frictionalPressure
        (
            const volScalarField& alphap,
            const volScalarField& alphaMax
        ) const;

        virtual tmp<volScalarField> frictionalPressurePrime
        (
            const volScalarField& alphap,
            const volScalarField& alphaMax
        ) const;

        virtual tmp<volScalarField> nu
        (
            const phaseModel& phase,
            const volScalarField& alphap,
            const volScalarField& alphaMax,
            const volScalarField& pf,
            const volSymmTensorField& D
        ) const;

        virtual tmp<volScalarField> alphaMinFriction
        (
            const volScalarField& alphap,
            const volScalarField& alphaMax
        ) const;

        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frictionalStressModels
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
