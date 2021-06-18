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
    Foam::diameterModels::qbmmDiameterModel

Description
    Diameter model using quadrature based moment methods

SourceFiles
    qbmmDiameterModel.C

\*---------------------------------------------------------------------------*/

#ifndef qbmmDiameterModel_H
#define qbmmDiameterModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "diameterModel.H"
#include "univariateODEPopulationBalance.H"

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                           Class qbmmDiameter Declaration
\*---------------------------------------------------------------------------*/

class qbmmDiameterModel
:
    public diameterModel
{
// Private data

    //- Population balance properties dictionary
    IOdictionary pbeDict_;

    //- QBMM model
    PDFTransportModels::populationBalanceModels::univariateODEPopulationBalance pbe_;

    //- Does the pbe use number density
    bool numberDensity_;

    //- Reference to volume fraction
    const volScalarField& vf_;

    //- Reference to volume fraction weighted diameter
    const volScalarField& vfD_;


public:

    //- Runtime type information
    TypeName("qbmm");

    // Constructor
    qbmmDiameterModel
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );


    //- Destructor
    virtual ~qbmmDiameterModel();


    // Member Functions

        //- Update fluxes
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Post update
        virtual void postUpdate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
