/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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
    Foam::multicomponentBlastThermo

Description
    Class to handle species transport

SourceFiles
    multicomponentBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentBlastThermo_H
#define multicomponentBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegrationSystem.H"
#include "PtrListDictionary.H"
#include "basicSpecieBlastMixture.H"
#include "fvScalarMatrix.H"
#include "speciesTable.H"
#include "fvModels.H"
#include "fvConstraints.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multicomponentBlastThermo Declaration
\*---------------------------------------------------------------------------*/

class multicomponentBlastThermo
:
    public basicSpecieBlastMixture
{
//- Private data

    const fvMesh& mesh_;

    //- Name of master thermo class
    word masterName_;

    //- Are the mass fractions normalized
    bool normalise_;

protected:
//- Protected data

    //- Mass transfer rates passed from an external source
    PtrListDictionary<volScalarField::Internal> massTransferRates_;

    //- Sources to solve implicitly
    PtrListDictionary<fvScalarMatrix> implicitSources_;

    //- Mass flux
    UautoPtr<const surfaceScalarField> alphaRhoPhiPtr_;

    //- Mass
    UautoPtr<const volScalarField> alphaRhoPtr_;

    class integrator
    :
        public timeIntegrationSystem
    {
    // Private data
        const fvMesh& mesh_;
        PtrList<volScalarField>& Y_;
        PtrListDictionary<volScalarField::Internal>& massTransferRates_;
        PtrListDictionary<fvScalarMatrix>& implicitSources_;
        const List<bool>& active_;
        const volScalarField& alphaRho_;
        const surfaceScalarField& alphaRhoPhi_;
        bool normalize_;

    public:
        integrator
        (
            const fvMesh& mesh,
            PtrList<volScalarField>& Y,
            PtrListDictionary<volScalarField::Internal>& massTransferRates,
            PtrListDictionary<fvScalarMatrix>& implicitSources,
            const List<bool>& active,
            const word& alphaRhoName,
            const word& alphaRhoPhiName,
            const bool normalize
        );

        virtual ~integrator();

        virtual void update();
        virtual void solve();
        virtual void postUpdate();

        bool normalize() const
        {
            return normalize_;
        }
    };

    autoPtr<integrator> integratorPtr_;

    // Return the mass field
    const volScalarField& alphaRho() const
    {
        return alphaRhoPtr_();
    }

    //- Return the mass flux
    const surfaceScalarField& alphaRhoPhi() const
    {
        return alphaRhoPhiPtr_();
    }


public:

    // Constructor
    multicomponentBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    // Constructor
    multicomponentBlastThermo
    (
        const speciesTable& species,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );


    //- Destructor
    virtual ~multicomponentBlastThermo();


    const fvMesh& mesh() const
    {
        return mesh_;
    }

    //- ODE functions

        //- Initialize models
        virtual void initializeModels();

        //- Pre update step
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Post update
        virtual void postUpdate();

        //- Add an external change to a thermodynamic field if present
        //  i.e. mass transfer
        virtual void clearDeltas();

        //- Add an external change to a thermodynamic field if present
        //  i.e. mass transfer
        virtual void clearSources();

        //- Add an external change to a thermodynamic field if present
        //  i.e. mass transfer
        virtual void addDelta(const word& name, tmp<volScalarField>&& delta);

        //- Add an external change to a thermodynamic field if present
        //  i.e. mass transfer
        virtual void addDelta(const word& name, tmp<volScalarField::Internal>&& delta);

        //- Add an external change to a thermodynamic field if present
        //  i.e. mass transfer
        virtual void addDelta(const word& name, const volScalarField::Internal& delta);

        //- Add an external change to a thermodynamic field if present
        //  i.e. mass transfer
        virtual void addSource(const word& name, tmp<fvScalarMatrix>& source);

        //- Add an external change to a thermodynamic field if present
        //  i.e. mass transfer
        virtual void addSource(const word& name, tmp<fvScalarMatrix>&& source);


     // Access

        //- Return if the composition is normalised
        bool normalise() const
        {
            return normalise_;
        }

        //- Access if the composition is normalised
        bool& normalise()
        {
            return normalise_;
        }

        //- Return the composition of the multi-component mixture
        virtual basicSpecieBlastMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicSpecieBlastMixture& composition() const = 0;

        //- Return the name of the thermo physics
        virtual word thermoName() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentBlastThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
