/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
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
    Foam::activationModel

Description
    Base model to describe activtion of an energetic material

SourceFiles
    activationModel.C
    newactivationModel.C

\*---------------------------------------------------------------------------*/

#ifndef activationModel_H
#define activationModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class activationModel Declaration
\*---------------------------------------------------------------------------*/

class activationModel
{
protected:

    // Protected data

        //- Progress variable
        volScalarField lambda_;

        //- Specific activation energy
        dimensionedScalar e0_;

        //- Exponent used to blend reacted and unreacted states (\lambda^m)
        scalar lambdaExp_;

        //- Density name
        word alphaRhoName_;

        //- Name of flux
        word alphaRhoPhiName_;

        //- Maximum change in lambda per time step
        scalar maxDLambda_;

        //- Is the change in lambda limited
        bool limit_;

        //- Stored old field indexes
        labelList oldIs_;
        label nOld_;

        //- Stored delta indexes
        labelList deltaIs_;
        label nDelta_;

        //- Stored old values
        PtrList<volScalarField> lambdaOld_;

        //- Stored changes in lambda
        tmp<volScalarField> ddtLambda_;

        //- Stored changes in lambda
        PtrList<volScalarField> deltaLambda_;

        //- Stored changes in lambda due to advection
        PtrList<volScalarField> deltaAlphaRhoLambda_;


    // Protected functions

        //- Limit change in lambda
        void limit();

        //- Return the time rate of chage of lambda
        virtual tmp<volScalarField> delta() const = 0;


public:

    //- Runtime type information
    ClassName("activationModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            activationModel,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const word& phaseName
            ),
            (mesh, dict, phaseName)
        );

    // Constructor
    activationModel
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );


    //- Destructor
    virtual ~activationModel();


    // Selectors

        static autoPtr<activationModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName
        );


    // Member Functions

        //- Solve sub-step stepi
        virtual void solve
        (
            const label stepi,
            const scalarList& ai,
            const scalarList& bi
        );

        //- Set old lists and fluxes (initialization of fields)
        virtual void setODEFields
        (
            const label nSteps,
            const labelList& oldIs,
            const label& nOld,
            const labelList& deltaIs,
            const label nDelta
        );

        //- Remove stored fields
        virtual void clearODEFields();

        //- Return the specific detonation energy
        const dimensionedScalar& e0() const
        {
            return e0_;
        }

        //- Return const reference to progress variable
        const volScalarField& lambda() const
        {
            return lambda_;
        }

        //- Return lambda to the m power for blending
        tmp<volScalarField> lambdaPow() const
        {
            return pow(lambda_, lambdaExp_);
        }

        //- Return ambda to the m power for blending for patchi
        tmp<scalarField> lambdaPow(const label patchi) const
        {
            return pow(lambda_.boundaryField()[patchi], lambdaExp_);
        }

        //- Return ambda to the m power for blending for celli
        scalar lambdaPowi(const label celli) const
        {
            return pow(lambda_[celli], lambdaExp_);
        }

        //- Return energy source
        virtual tmp<volScalarField> ddtLambda() const;

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
