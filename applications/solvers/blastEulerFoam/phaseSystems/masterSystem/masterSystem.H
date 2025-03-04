/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::masterSystem

Description
    Class to handle a group of systems that need to be updated once per group
SourceFiles
    masterSystem.C

\*---------------------------------------------------------------------------*/

#ifndef masterSystem_H
#define masterSystem_H

#include "phaseSystem.H"
#include "phaseModel.H"
#include "UPtrList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class packingLimitModel;

/*---------------------------------------------------------------------------*\
                 Class kineticTheorySystem Declaration
\*---------------------------------------------------------------------------*/

class masterSystem
:
    public regIOobject
{
protected:
    // Protected data

        //- Constructing dictionary
        const dictionary& dict_;

        //- Group of the system
        word group_;

        //- Write totals
        bool writeTotal_;

        //- Const refrence to multiphaseSystem
        const phaseSystem& fluid_;

        //- List of phases part of the system
        UPtrList<phaseModel> phases_;

        //- List of granular phase indexes
        labelList phaseIndexes_;


        // Total fields

            //- Total phase-fraction
            autoPtr<volScalarField> alphaPtr_;

            //- Phase averaged density
            mutable autoPtr<volScalarField> rhoPtr_;

            //- Phase averaged velocity
            autoPtr<volVectorField> UPtr_;

            //- Phase averaged volumetric flux
            mutable autoPtr<surfaceScalarField> phiPtr_;

            //- Phase averaged volume fraction flux
            mutable autoPtr<surfaceScalarField> alphaPhiPtr_;

            //- Was residual alpha read
            bool readResidualAlpha_;

            //- Residual phase fraction
            dimensionedScalar residualAlpha_;

            //- Was residual alpha read
            bool readResidualRho_;

            //- Residual phase density
            dimensionedScalar residualRho_;

        // Packing limits (optional)

            //- Run-time selected packing limit model
            autoPtr<packingLimitModel> packingLimitModel_;

            //- Alpha max field
            volScalarField alphaMax_;

            //- Minimum phase alpha max
            scalar minAlphaMax_;


public:

    TypeName("masterSystem");

    // Constructor
    masterSystem
    (
        const word& type,
        const phaseSystem&,
        const dictionary&
    );

    //- Lookup or construct kineticTheorySystem
    template<class System>
    static System& New(const phaseSystem&);

    //- Return a new object pointer
    template<class System>
    static autoPtr<System> NewPtr(const phaseSystem&);

    //- Destructor
    ~masterSystem();

    //- Disallow default bitwise copy construct
    masterSystem(const masterSystem&) = delete;


    // Member Functions

        //- Return the group name
        const word& group() const
        {
            return group_;
        }

        //- Return fluid
        const phaseSystem& fluid() const
        {
            return fluid_;
        }

        //- Return if the phase is polydisperse
        bool polydisperse() const
        {
            return phases_.size() > 1;
        }

        //- Return the list of phases
        const UPtrList<phaseModel>& phases() const
        {
            return phases_;
        }

        //- Const refrence to list of granular indexes
        const labelList& phaseIndexes() const;

        //- Add phase to list of kinetic theory phases.
        virtual void addPhase(phaseModel&);

        //- Return true if phase is in the system
        bool contains(const phaseModel&) const;

        //- Return true if phase name is in the system
        bool contains(const word&) const;

        //- Return true if phase is in the system
        label whichPhase(const phaseModel&) const;

        //- Return true if phase name is in the system
        label whichPhase(const word&) const;

        //- Return residual alpha
        const dimensionedScalar& residualAlpha() const
        {
            return residualAlpha_;
        }

        //- Return residual rho
        const dimensionedScalar& residualRho() const
        {
            return residualRho_;
        }

        //- Return total volume fraction
        const volScalarField& alpha() const;

        //- Return phase averaged density
        const volScalarField& rho() const;

        //- Return phase averaged velocity
        const volVectorField& U() const;

        //- Return phase averaged volumetric flux
        const surfaceScalarField& phi() const;

        //- Return total volume fraction flux
        tmp<surfaceScalarField> alphaPhi() const;

        //- Return the minimum max alpha
        scalar minAlphaMax() const
        {
            return minAlphaMax_;
        }

        //- Return maximum volume fraction
        virtual const volScalarField& alphaMax() const
        {
            return alphaMax_;
        }

        //- Return maximum volume fraction
        virtual scalar alphaMax(const label celli) const
        {
            return alphaMax_[celli];
        }

        //- Initialize the model after everything is constructed
        virtual void initialize()
        {}

        //- Correct total volume fraction
        virtual void correctAlpha();

        //- Update
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve()
        {}

        //- Post update
        virtual void postUpdate()
        {}

        //- Remove stored fields
        virtual void clear()
        {}

        //- Write using setting from DB
        virtual bool write(const bool write = true) const
        {
            return true;
        }

        //- Write using given format, version and compression
        virtual bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType,
            const bool write
        ) const
        {
            return true;
        }

        virtual bool writeData(Ostream&) const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const masterSystem&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "masterSystemTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
