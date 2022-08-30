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


protected:
//- Protected data

    //- Mass transfer rates passed from an external source
    PtrListDictionary<volScalarField> massTransferRates_;

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
        PtrListDictionary<volScalarField>& massTransferRates_;
        PtrListDictionary<fvScalarMatrix>& implicitSources_;
        const List<bool>& active_;
        const volScalarField& alphaRho_;
        const surfaceScalarField& alphaRhoPhi_;

    public:
        integrator
        (
            const fvMesh& mesh,
            PtrList<volScalarField>& Y,
            PtrListDictionary<volScalarField>& massTransferRates,
            PtrListDictionary<fvScalarMatrix>& implicitSources,
            const List<bool>& active,
            const word& alphaRhoName,
            const word& alphaRhoPhiName
        );

        virtual ~integrator();

        virtual void update();
        virtual void solve();
        virtual void postUpdate();
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

    //- Correct mass fractions
    void correct();



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

        //- Update the mixtures
        virtual void updateMixture() = 0;

        //- Pre update step
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Post update
        virtual void postUpdate();

        //- Add an external change to a thermodynamic field if present
        //  i.e. mass transfer
        virtual void addDelta(const word& name, tmp<volScalarField>& delta);

        //- Add an external change to a thermodynamic field if present
        //  i.e. mass transfer
        virtual void addDelta(const word& name, const volScalarField::Internal& delta);

        //- Add an external change to a thermodynamic field if present
        //  i.e. mass transfer
        virtual void addSource(const word& name, tmp<fvScalarMatrix>& source);


     // Access

        //- Return the composition of the multi-component mixture
        virtual basicSpecieBlastMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicSpecieBlastMixture& composition() const = 0;

        //- Flame temperature [K]
        virtual scalar flameT(const label speciei) const = 0;

        //- Return the pressure of a cell
        virtual scalar cellp(const label celli) const = 0;

        //- Pressure [kg/m/s^2]
        virtual scalar p
        (
            const label speciei,
            const scalar rho,
            const scalar e,
            const scalar T
        ) const = 0;

        //- Density [kg/m^3]
        virtual tmp<volScalarField> p
        (
            const label speciei,
            const volScalarField& rho,
            const volScalarField& e,
            const volScalarField& T
        ) const = 0;

        //- Return derivative of pressure w.r.t. density
        virtual scalar dpdRho
        (
            const label speciei,
            const scalar rho,
            const scalar e,
            const scalar T
        ) const = 0;

        //- Return derivative of pressure w.r.t. Temperature
        virtual scalar dpdT
        (
            const label speciei,
            const scalar rho,
            const scalar e,
            const scalar T
        ) const = 0;

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
