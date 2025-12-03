/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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
    Foam::embeddedTimeIntegrator

Description
    Base class for time integration with embedded error

SourceFiles
    embeddedTimeIntegrator.C
    embeddedTimeIntegratorNew.C

\*---------------------------------------------------------------------------*/

#ifndef embeddedTimeIntegrator_H
#define embeddedTimeIntegrator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegrator.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class embeddedTimeIntegrator Declaration
\*---------------------------------------------------------------------------*/

class embeddedTimeIntegrator
:
    public timeIntegrator
{
protected:
// Protected variables

    // Step size controls
    dictionary dict_;
    scalar safeScale_;
    scalar alphaInc_;
    scalar alphaDec_;
    scalar minScale_;
    scalar maxScale_;

    //- New time step size
    scalar deltaT_;

    //- Adjust deltaT
    bool adjust_;



    // Protected member functions

        //- Initialize ODE sizes
        virtual void initialize();


        // Constructor
        embeddedTimeIntegrator
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const bool embedded
        );

public:

    //- Runtime type information
    TypeName("embeddedTimeIntegrator");

    // Constructor
    embeddedTimeIntegrator
    (
        const objectRegistry& obr,
        const dictionary& dict
    );


    //- Destructor
    virtual ~embeddedTimeIntegrator();


    // Member Functions

        // Access

            //- Return new time step size
            scalar deltaT() const
            {
                return deltaT_;
            }

            //- Access if time step adjustment is allowed
            bool& adjust()
            {
                return adjust_;
            }


        // Edit

            //- Add a system to the
            virtual void addSystem(timeIntegrationSystemBase& system);

            //- Integrate with adaptive time stepping
            using timeIntegrator::integrate;
            virtual void integrate
            (
                const bool doExplicit,
                const bool doSave,
                const bool doImplicit,
                const bool doPost,
                const bool doClear
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
