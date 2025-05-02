/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
09-09-2024 Synthetik Applied Technologies: |    Calculate Courant number with
                                                blastFoam thermo
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::functionObjects::blastCourantNo

Description
    Calculates and writes the Mach number as a volScalarField.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    blastCourantNo.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_blastCourantNo_H
#define functionObjects_blastCourantNo_H

#include "fieldExpression.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class MachNo Declaration
\*---------------------------------------------------------------------------*/

class blastCourantNo
:
    public fieldExpression
{
    // Private Member data

        //- Name of phase
        word phaseName_;

        //- Name of velocity field
        word UName_;


    // Protected member functions

        virtual bool calc();


public:

    //- Runtime type information
    TypeName("blastCourantNo");

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        blastCourantNo
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        blastCourantNo(const blastCourantNo&) = delete;


    //- Destructor
    virtual ~blastCourantNo();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const blastCourantNo&) = delete;

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return {UName_};
        }
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
