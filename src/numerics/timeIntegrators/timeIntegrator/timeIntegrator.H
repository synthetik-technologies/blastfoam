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
    Foam::timeIntegrator

Description
    Base class for time integration

SourceFiles
    timeIntegrator.C
    newTimeIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef timeIntegrator_H
#define timeIntegrator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "runTimeSelectionTables.H"
#include "timeIntegratorCoeffs.H"
#include "Time.H"
#include "PtrList.H"
#include "scalarField.H"
#include "scalarList.H"

namespace Foam
{

class timeIntegrationSystemBase;

/*---------------------------------------------------------------------------*\
                           Class timeIntegrator Declaration
\*---------------------------------------------------------------------------*/

class timeIntegrator
:
    public regIOobject
{
protected:
// Protected variables

    //- Reference to mesh
    const objectRegistry& obr_;

    //- Ode coefficients
    autoPtr<timeIntegratorCoeffs> coeffs_;

    //- Reference to compressible system
    UPtrList<timeIntegrationSystemBase> systems_;

    //- Number of steps
    label nSteps_;

    //- Current step index
    label stepi_;

    //- Old value coefficients
    List<List<scalar>> as_;

    //- Delta coefficients
    List<List<scalar>> bs_;

    //- Time step fractions
    scalarList f_;

    //- Starting time step fractions
    scalarList f0_;

    //- Stored old field indexes
    labelList oldIs_;

    //- Number of stored fields
    label nOld_;

    //- Stored delta indexes
    labelList deltaIs_;

    //- Number of stored deltas
    label nDelta_;

    //- Current time index
    label curTimeIndex_;

    //- Is this a restarted time step
    bool restart_;


    // Protected member functions

        // Update coefficients
        virtual void update();

        //- Initialize ODE sizes
        virtual void initialize();


        // Constructor
        timeIntegrator
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const bool embedded
        );

public:

    //- Runtime type information
    TypeName("timeIntegrator");

    // Constructor
    timeIntegrator(const objectRegistry& obr, const dictionary& dict);


    //- Destructor
    virtual ~timeIntegrator();


    // Member Functions

        // Access

            //- Return the Time object
            inline const Time& time() const;

            //- Return if on a restarted time step
            inline bool restart() const;

            //- Return old coefficients
            inline const List<scalarList>& as() const;

            //- Return delta coefficients
            inline const List<scalarList>& bs() const;

            //- Return old coefficients for the current step
            inline const scalarList& a() const;

            //- Return delta coefficients for the current step
            inline const scalarList& b() const;

            //- Return the number of steps
            inline label nSteps() const;

            //- Return number of stored fields
            inline label nOld() const;

            //- Return the index in the list of old fields
            inline label getOldIndex(const label stepi) const;

            //- Return number of stored deltas
            inline label nDelta() const;

            //- Return the index in the list of delta fields
            inline label getDeltaIndex(const label stepi) const;

            //- Return current step
            inline label step() const;

            //- Return local step fraction
            inline scalar f() const;

            //- Return start fraction of the current step
            inline scalar f0() const;

            //- Is the current step the last
            inline bool firstStep() const;

            //- Is the current step the last
            inline bool finalStep() const;

            //- Return map for saving old values
            inline labelList oldSaveMap() const;

            //- Return map for saving delta values
            inline labelList deltaSaveMap() const;

            //- Return max Courant number
            inline scalar maxCo() const;


        // Edit

            //- Add a system to the
            virtual void addSystem(timeIntegrationSystemBase& system);

            //- Update all stored systems
            virtual void updateAll();

            //- Post update
            virtual void preUpdate();

            //- Integrate fluxes in time
            virtual void integrate();
            virtual void integrate
            (
                const bool doExplicit,
                const bool doSave,
                const bool doImplicit,
                const bool doPost,
                const bool doClear
            );

            //- Save explicit changes
            virtual void storeExplicit();

            //- Solve implicit
            virtual void solveExplicit();

            //- Solve implicit
            virtual void solveImplicit();

            //- Post update
            virtual void postUpdate();

            //- Clear all ODE fields (old/delta)
            virtual void clear();

            //- Reset fields to the old time
            virtual void reset();

            //- Return current step old volume
            virtual tmp<scalarField> V0() const;

            //- Return current step volume
            virtual tmp<scalarField> V() const;

            //- Total old volume
            virtual scalar totalV0() const;

            //- Total volume
            virtual scalar totalV() const;

        //- Dummy IO operations
        virtual bool writeData(Ostream& os) const
        {
            return true;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegratorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
