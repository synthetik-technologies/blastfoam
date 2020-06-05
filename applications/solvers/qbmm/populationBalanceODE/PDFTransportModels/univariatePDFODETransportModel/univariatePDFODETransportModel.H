/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 Alberto Passalacqua
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
    Foam::univariatePDFODETransportModel

Description
    Solve a univariate PDF transport equation with the extended quadrature
    method of moments.

SourceFiles
    univariatePDFODETransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef univariatePDFODETransportModel_H
#define univariatePDFODETransportModel_H

#include "PDFODETransportModel.H"
#include "quadratureApproximations.H"
#include "univariateMomentAdvection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{

/*---------------------------------------------------------------------------*\
                Class univariatePDFODETransportModel Declaration
\*---------------------------------------------------------------------------*/

class univariatePDFODETransportModel
:
    public PDFODETransportModel
{
protected:

    // Protected data

        //- Univariate quadrature approximation used to solve PDF transport
        scalarQuadratureApproximation quadrature_;

        //- Moment advection
        autoPtr<univariateMomentAdvection> momentAdvection_;

    // Protected member functions

        //- Return number of moments
        virtual label nMoments() const
        {
            return quadrature_.nMoments();
        }

        //- Calculate implicit source terms
        virtual tmp<fvScalarMatrix> implicitMomentSource
        (
            const volScalarMoment& moment
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
        univariatePDFODETransportModel
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh,
            const surfaceScalarField& phi,
            const word& support
        );


    //- Destructor
    virtual ~univariatePDFODETransportModel();


    // Member Functions

        //- Access

            //- Returns a const reference to the quadrature approximation
            inline const scalarQuadratureApproximation& quadrature() const;

            //- Return the maximum Courant number ensuring moment realizability
            inline scalar realizableCo() const;

        //-Edit

            //- Solve sub-step stepi
            virtual void solve
            (
                const label stepi,
                const scalarList& ai,
                const scalarList& bi
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PDFTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariatePDFODETransportModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
