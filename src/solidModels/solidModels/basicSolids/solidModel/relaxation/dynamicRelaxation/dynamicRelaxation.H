/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    dynamicRelaxation

Description
    Run-time selectable stabilisation (diffusion) term for the solid mechanics
    momentum equation.

Author
    Philip Cardiff, UCD.

SourceFiles
    dynamicRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicRelaxation_H
#define dynamicRelaxation_H

#include "dictionary.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dynamicRelaxation Declaration
\*---------------------------------------------------------------------------*/

class dynamicRelaxation
{
private:
    // Private data

        //- Is relaxation active
        bool relax_;

        //- Start time of relaxation
        scalar start_;

        //- Number of steps before scaling
        label nSteps_;

        //- Current time index
        label curIndex_;

        //- Stored kinetic energies
        DynamicList<scalar> Ks_;

        //- Scaling factor
        scalar scale_;

        //- Optional probe location
        bool useProbe_;
        vector probe_;

        //- Disallow default bitwise copy construct
        dynamicRelaxation(const dynamicRelaxation&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicRelaxation&);

public:

    //- Runtime type information
    TypeName("dynamicRelaxation");


    // Constructors

        //- Construct from dictionary
        dynamicRelaxation(const dictionary& dict);


    // Destructor

        virtual ~dynamicRelaxation()
        {}


    // Member Functions

        //- Read parameters
        void read(const dictionary&);

        //- Relax U field
        void relax(volVectorField& U, const volScalarField& rho);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
