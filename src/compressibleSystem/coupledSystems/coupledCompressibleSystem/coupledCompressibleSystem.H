/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
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
    Foam::coupledCompressibleSystem

Description
    Compressible system containing more than two equation of states.

SourceFiles
    coupledCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef coupledCompressibleSystem_H
#define coupledCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "compressibleSystem.H"
#include "UautoPtr.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class coupledCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class coupledCompressibleSystem
{
protected:
    // Protected data

        //- Total volume fraction
        volScalarField volumeFraction_;

        //- Pointer to dispersed phase volume fraction
        UautoPtr<const volScalarField> alphadPtr_;

        //- Phasic density
        volScalarField alphaRho_;

public:

    TypeName("coupledCompressibleSystem");

    // Constructor
    coupledCompressibleSystem
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    // Selector
    static autoPtr<compressibleSystem> New
    (
        const fvMesh& mesh
    );

    //- Destructor
    virtual ~coupledCompressibleSystem();


    // Member functions

        //- Set a pointer to the dispersed phase volume fraction
        virtual void setDispersedVolumeFraction(const volScalarField& alphad)
        {
            alphadPtr_.set(&alphad);
        }

        //- Return const access to volume fraction
        virtual const volScalarField& volumeFraction() const
        {
            return volumeFraction_;
        }

        //- Non-cont access to volume fraction
        virtual volScalarField& volumeFraction()
        {
            return volumeFraction_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
