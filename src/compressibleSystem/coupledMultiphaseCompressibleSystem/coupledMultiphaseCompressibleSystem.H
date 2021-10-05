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

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class coupledMultiphaseCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class coupledMultiphaseCompressibleSystem
:
    public multiphaseCompressibleSystem
{
    // Private data

        //- Total volume fraction
        volScalarField volumeFraction_;

        //- Pointer to dispersed phase volume fraction
        UautoPtr<const volScalarField> alphadPtr_;

        //- Phasic density
        volScalarField alphaRho_;


public:

    TypeName("coupledMultiphaseCompressibleSystem");

    // Constructor
    coupledMultiphaseCompressibleSystem(const fvMesh& mesh);


    //- Destructor
    ~coupledMultiphaseCompressibleSystem();


    // Member functions

        //- Set a pointer to the dispersed phase volume fraction
        void setDispersedVolumeFraction(const volScalarField& alphad)
        {
            alphadPtr_.set(&alphad);
            encode();
        }

        //- Return const access to volume fraction
        const volScalarField& volumeFraction() const
        {
            return volumeFraction_;
        }

        //- Non-cont access to volume fraction
        volScalarField& volumeFraction()
        {
            return volumeFraction_;
        }

        //- Decode primative variables
        virtual void decode();

        //- Update fluxes
        virtual void update();

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
