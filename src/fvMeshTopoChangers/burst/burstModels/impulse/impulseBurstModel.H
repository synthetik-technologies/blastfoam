/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022-2025
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
    Foam::burstModels::impulse

Description
    Base class for time integration

SourceFiles
    impulseBurstModel.C
    impulseBurstModel.C

\*---------------------------------------------------------------------------*/

#ifndef impulseBurstModel_H
#define impulseBurstModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fieldBurstModel.H"


namespace Foam
{
namespace burstModels
{

/*---------------------------------------------------------------------------*\
                           Class burstModel Declaration
\*---------------------------------------------------------------------------*/

class impulse
:
    public field
{
public:

    //- Runtime type information
    TypeName("impulse");

    //- Constructor
    impulse
    (
        const dictionary& dict,
        const bool coupled
    );

    //- Clone
    virtual autoPtr<burstModel> clone() const
    {
        return autoPtr<burstModel>(new impulse(*this));
    }


    //- Destructor
    virtual ~impulse();


    // Member Functions

        //- Write
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace burstModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
