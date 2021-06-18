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
    Foam::diameterModels::reactingDiameterModel

Description
    Reacting diameter model used for mass transfer

SourceFiles
    reactingDiameterModel.C

\*---------------------------------------------------------------------------*/

#ifndef reactingDiameterModel_H
#define reactingDiameterModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "diameterModel.H"
#include "reactionRate.H"

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                           Class reactingDiameterModel Declaration
\*---------------------------------------------------------------------------*/

class reactingDiameterModel
:
    public diameterModel
{
// Private data

    //- Reaction rate model
    autoPtr<reactionRate> rate_;

    //- Name of pressure field
    const word pName_;

    //- Name of temperature field
    const word TName_;

    //- Stored change in volume per unit time
    volScalarField dVdt_;

    //- Does the diameter field need to be read
    virtual bool requireRead()
    {
        return true;
    }

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

    //- Return if the pressure model requires inputs
    virtual bool requireInput() const
    {
        return true;
    }

    // ODE functions

        //- Solve sub-step stepi
        virtual void solve();

        //- Solve sub-step stepi given pressure and temperature fields
        virtual void solve(const volScalarField& p, const volScalarField& T);

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
