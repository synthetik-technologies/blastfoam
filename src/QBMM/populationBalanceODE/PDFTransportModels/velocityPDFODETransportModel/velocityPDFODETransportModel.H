/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 Alberto Passalacqua
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
    Foam::velocityPDFODETransportModel

Description
    Solve a velocity PDF transport equation with the hyperbolic conditional
    quadrature method of moments

SourceFiles
    velocityPDFODETransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef velocityPDFODETransportModel_H
#define velocityPDFODETransportModel_H

#include "PDFODETransportModel.H"
#include "quadratureApproximations.H"
#include "velocityMomentAdvection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{

/*---------------------------------------------------------------------------*\
                Class velocityPDFODETransportModel Declaration
\*---------------------------------------------------------------------------*/

class velocityPDFODETransportModel
:
    public PDFODETransportModel
{
protected:

    // Protected data

        //- velocity quadrature approximation used to solve PDF transport
        velocityQuadratureApproximation quadrature_;

        //- Moment advection
        autoPtr<velocityMomentAdvection> momentAdvection_;

    // Protected member functions

        //- Return number of moments
        virtual label nMoments() const
        {
            return quadrature_.nMoments();
        }

        //- Update implicit collison source terms
        virtual void updateImplicitMomentSource() = 0;

        //- Calculate implicit source terms
        virtual tmp<fvScalarMatrix> implicitMomentSource
        (
            const volVelocityMoment& moment
        ) = 0;

        //- Calculate explicit source terms with realizable ODE solver
        virtual void explicitMomentSource() = 0;

        //- Are moments updated due to sources
        virtual bool solveMomentSources() const = 0;

        //- Is the realizable ode solver used to solve sources
        virtual bool solveMomentOde() const = 0;


public:

    // Constructors

        //- Construct from components
        velocityPDFODETransportModel
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh,
            const word& support
        );


    //- Destructor
    virtual ~velocityPDFODETransportModel();


    // Member Functions

        //- Access

            //- Returns a const reference to the quadrature approximation
            virtual inline const velocityQuadratureApproximation& quadrature() const;

            //- Return the maximum Courant number ensuring moment realizability
            virtual inline scalar realizableCo() const;

            //- Return the maximum Courant number based on velocities
            virtual inline scalar CoNum() const;


        //-Edit

            //- Update fluxes
            virtual void update();

            //- Solve sub-step stepi
            virtual void solve
            (
                const label stepi,
                const scalarList& ai,
                const scalarList& bi
            );

            //- Solve implicit terms and ode
            virtual void postUpdate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PDFTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "velocityPDFODETransportModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
