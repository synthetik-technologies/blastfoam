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
    Foam::coupledMultiphaseCompressibleSystem

Description
    Compressible system containing more than two equation of states.

SourceFiles
    coupledMultiphaseCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef coupledMultiphaseCompressibleSystem_H
#define coupledMultiphaseCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiphaseCompressibleSystem.H"
#include "coupledCompressibleSystem.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class coupledMultiphaseCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class coupledMultiphaseCompressibleSystem
:
    public coupledCompressibleSystem,
    public multiphaseCompressibleSystem
{
    // Private data

        //- Calculate new alpha and rho fields
        virtual void calcAlphas();

public:

    TypeName("multiphaseCoupled");

    // Constructor
    coupledMultiphaseCompressibleSystem(const fvMesh& mesh);


    //- Destructor
    ~coupledMultiphaseCompressibleSystem();


    // Member functions

        //- Set a pointer to the dispersed phase volume fraction
        virtual void setDispersedVolumeFraction(const volScalarField& alphad)
        {
            coupledCompressibleSystem::setDispersedVolumeFraction(alphad);
            encode();
        }

        //- Return const access to volume fraction
        virtual const volScalarField& volumeFraction() const
        {
            return coupledCompressibleSystem::volumeFraction();
        }

        //- Non-cont access to volume fraction
        virtual volScalarField& volumeFraction()
        {
            return coupledCompressibleSystem::volumeFraction();
        }

        //- Const access to total mass
        virtual const volScalarField& rhoEff() const
        {
            return alphaRho_;
        }

        //- Access to total mass
        virtual volScalarField& rhoEff()
        {
            return alphaRho_;
        }

        //- Decode primative variables
        virtual void decode();

        //- Encode conserved variables
        virtual void encode();

        //- Solve a sub-step
        virtual void solve();

        //- Solve implicit terms
        virtual void postUpdate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
