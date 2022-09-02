/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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
    Foam::burstModels::pressureImpulse

Description
    Burst model using pressure and impulse

SourceFiles
    pressureImpulseBurstModel.C

\*---------------------------------------------------------------------------*/

#ifndef overpressureImpulseBurstModel_H
#define overpressureImpulseBurstModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pressureBurstModel.H"


namespace Foam
{
namespace burstModels
{

/*---------------------------------------------------------------------------*\
                    Class pressureImpulse Declaration
\*---------------------------------------------------------------------------*/

class pressureImpulse
:
    public pressure
{
protected:
// Protected variables

    //- Name of impulse field
    word impulseName_;

    //- Impulse that the patch will burst
    scalar impulseBurst_;


public:

    //- Runtime type information
    TypeName("pressureImpulse");

    // Constructor
    pressureImpulse(const dictionary& dict);

    //- Clone
    virtual autoPtr<burstModel> clone() const
    {
        return autoPtr<burstModel>(new pressureImpulse(*this));
    }

    //- Destructor
    virtual ~pressureImpulse();


    // Member Functions

        //- Return if partial burst is used
        inline const word& impulseName() const
        {
            return impulseName_;
        }


        //- Update the given intact field
        virtual bool update
        (
            const fvPatch&,
            scalarField& intact
        ) const;

        //- Write
        virtual void writeData(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace burstModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
