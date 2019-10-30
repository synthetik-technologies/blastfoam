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
    Foam::equationOfState

Description

SourceFiles
    equationOfState.C
    newEquationOfState.C

\*---------------------------------------------------------------------------*/

#ifndef equationOfState_H
#define equationOfState_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class equationOfState Declaration
\*---------------------------------------------------------------------------*/

class equationOfState
:
    public volScalarField
{
protected:

    // Protected data

        //- Dictionary containing coefficients
        const dictionary& eosDict_;

        //- Constant reference o internal energy
        const volScalarField& e_;

        //- Phase density
        volScalarField rho_;

        //- Volume fraction weighted density
        volScalarField alphaRho_;

        //- Volume fraction flux
        surfaceScalarField alphaPhi_;

        //- Mass flux
        surfaceScalarField alphaRhoPhi_;

        //- Minimum volume fraction
        dimensionedScalar residualAlpha_;

        //- Minimum density
        dimensionedScalar residualRho_;


public:

    //- Runtime type information
    ClassName("equationOfState");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            equationOfState,
            dictionary,
            (
                const volScalarField& e,
                const dictionary& dict
            ),
            (e, dict)
        );

    // Constructor
    equationOfState(const volScalarField& e, const dictionary& dict);


    //- Destructor
    virtual ~equationOfState();


    // Selectors

        static autoPtr<equationOfState> New
        (
            const volScalarField& e,
            const dictionary& dict
        );


    // Member Functions

        //- Update
        virtual void update()
        {}

        //- Return pressure
        virtual tmp<volScalarField> pByGamma() const;

        //- Return density with minimum value
        tmp<volScalarField> stabilizeRho() const
        {
            return Foam::max(rho_, residualRho_);
        }

        //- Reference internal energy
        virtual tmp<volScalarField> e(const volScalarField& p) const;

        //- Gamma function
        virtual tmp<volScalarField> Gamma() const = 0;


        //- Pi function
        virtual tmp<volScalarField> Pi() const = 0;

        //- Delta coefficient
        virtual tmp<volScalarField> delta() const = 0;

        //- Scaling coefficient
        tmp<volScalarField> xi() const;

        //- Speed of sound squared
        tmp<volScalarField> cSqr(const volScalarField& P) const;

        // Access functions

            //- Return residual volume fraction
            const dimensionedScalar& residualAlpha() const
            {
                return residualAlpha_;
            }

            //- Return residual density
            const dimensionedScalar& residualRho() const
            {
                return residualRho_;
            }

            //- Return constant reference to density field
            const volScalarField& rho() const
            {
                return rho_;
            }

            //- Return non-constant reference to density field
            volScalarField& rho()
            {
                return rho_;
            }

            //- Return constant reference to mass
            const volScalarField& alphaRho() const
            {
                return alphaRho_;
            }

            //- Return non-constant reference to mass
            volScalarField& alphaRho()
            {
                return alphaRho_;
            }

            //- Return constant reference to volume fraction flux
            const surfaceScalarField& alphaPhi() const
            {
                return alphaPhi_;
            }

            //- Return non-constant reference to volume fraction flux
            surfaceScalarField& alphaPhi()
            {
                return alphaPhi_;
            }

            //- Return constant reference to mass flux
            const surfaceScalarField& alphaRhoPhi() const
            {
                return alphaRhoPhi_;
            }

            //- Return non-constant reference to mass flux
            surfaceScalarField& alphaRhoPhi()
            {
                return alphaRhoPhi_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
