/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technologies
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
    Foam::activationModels::programmedIgnitionActivation

Description
    Activation model using linear velocity and compression based ignition options

    References:
    \verbatim
        Giroux, E. D. (1973).
        Hemp User's Manual,
        Lawrence Livermore Laboratory
    \endverbatim

SourceFiles
    programmedIgnitionActivation.C

\*---------------------------------------------------------------------------*/

#ifndef programmedIgnitionActivation_H
#define programmedIgnitionActivation_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "activationModel.H"

namespace Foam
{
namespace activationModels
{

/*---------------------------------------------------------------------------*\
                           Class programmedIgnitionActivation Declaration
\*---------------------------------------------------------------------------*/

class programmedIgnitionActivation
:
    public activationModel
{
public:
    //- Burn models
        enum burnModel
        {
            BETA,           // Compression based
            PROGRAMMED,     // Programmed burn
            PROGRAMMEDBETA  // Compression and programmed burn
        };

        static const NamedEnum<burnModel, 3> burnModelNames_;

private:
// Private data

        //- Const reference to density
        const volScalarField& rho_;

        //- Detonation velocity
        dimensionedScalar vDet_;

        //- Chapmann-Jouguet pressure
        dimensionedScalar Pcj_;

        //- Reference density
        dimensionedScalar rho0_;

        //- Chapmann-Jouguet specific volume
        dimensionedScalar Vcj_;

        //- Burn model
        burnModel model_;

        //- Ignition time
        volScalarField::Internal tIgn_;

        //- Return the time rate of chage of lambda
        virtual tmp<volScalarField> delta() const;

        //- Correct lambda field after advection
        virtual void correct();

        //- Do detonation points need to be specified
        virtual bool needDetonationPoints() const
        {
            return false;
        }



public:

    //- Runtime type information
    TypeName("programmedIgnition");

    // Constructor
    programmedIgnitionActivation
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );


    //- Destructor
    virtual ~programmedIgnitionActivation();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace activationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
