/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020
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
    Foam::diameterModels::reactingDiameterModel

Description
    Reacting diameter model used for mass transfer

    With an optional oxidiser phase the burn rate is scaled by the fraction
    of the surrounding fluid volume that this phase occupies and, with
    oxidiser species, by min(1, sum of their mass fractions/YoxidiserRef):

    \verbatim
        oxidiserPhase   gas;
        oxidiserSpecies (O2 CO2 H2O);   // optional
        YoxidiserRef    0.24;           // full rate at this mass fraction
    \endverbatim

    Without oxidiserPhase the particles burn wherever the rate model allows.

SourceFiles
    reactingDiameterModel.C

\*---------------------------------------------------------------------------*/

#ifndef reactingDiameterModel_H
#define reactingDiameterModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "diameterModel.H"
#include "diameterReactionRate.H"

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                           Class reactingDiameterModel Declaration
\*---------------------------------------------------------------------------*/

class reactingDiameterModel
:
    public diameter
{
// Private data

    //- Reaction rate model
    autoPtr<diameterReactionRate> rate_;

    //- Name of pressure field
    const word pName_;

    //- Name of temperature field
    const word TName_;

    //- Phase carrying the oxidiser (empty: no oxidiser limit)
    const word oxidiserPhase_;

    //- Oxidiser species in the oxidiser phase (empty: the whole phase)
    const wordList oxidiserSpecies_;

    //- Oxidiser mass fraction that gives the full burn rate
    const scalar YoxidiserRef_;

    //- Number density
    volScalarField N_;

    //- Initial diameter field
    volScalarField d0_;

    //- Stored change in volume per unit time
    volScalarField dVdt_;


    // Private Member Functions

        //- Burn rate scaling by the available oxidiser, in [0, 1]
        tmp<volScalarField> oxidiserFactor() const;

public:

    //- Runtime type information
    TypeName("reacting");

    // Constructor
    reactingDiameterModel
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );


    //- Destructor
    virtual ~reactingDiameterModel();

    // ODE functions

        //- Update number density
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Solve sub-step stepi given pressure and temperature fields
        virtual void solve(const volScalarField& p, const volScalarField& T);

        //- Return initial unburned diameter field
        virtual const volScalarField& d0() const
        {
            return d0_;
        }

        //- Return initial unburned diameter for a cell
        virtual scalar celld0(const label celli) const
        {
            return d0_[celli];
        }

        //- Return change in radius per unit time
        //  Only used in mass transfer
        virtual tmp<volScalarField> dMdt() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
