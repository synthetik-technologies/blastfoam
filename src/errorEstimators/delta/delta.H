/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2020
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
    Foam::errorEstimators::delta

Description
    Error is defined by the relative difference of a field over a face

SourceFiles
    delta.C

\*---------------------------------------------------------------------------*/

#ifndef delta_H
#define delta_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "errorEstimator.H"

namespace Foam
{
namespace errorEstimators
{

/*---------------------------------------------------------------------------*\
                           Class delta Declaration
\*---------------------------------------------------------------------------*/

class delta
:
    public errorEstimator
{
// Private data

        //- Field name
        const word fieldName_;

        //- Minimum value to scale by
        scalar minVal_;

        //- Stored pointer to error function
        scalar (*errorFunc_)(const scalar&, const scalar&, const scalar&);


        // Functions to calculate error

            //- Unscaled error
            static scalar unscaledError
            (
                const scalar& own,
                const scalar& nei,
                const scalar& minVal
            )
            {
                return mag(own - nei);
            }

            //- Scaled error
            static scalar scaledError
            (
                const scalar& own,
                const scalar& nei,
                const scalar& minVal
            )
            {
                return mag(own - nei)/max(min(mag(own), mag(nei)), minVal);
            }


public:

    //- Runtime type information
    TypeName("delta");

    // Constructor
    delta
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& name
    );


    //- Destructor
    virtual ~delta();


    // Member Functions

        //- Gamma function
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace errorEstimators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
