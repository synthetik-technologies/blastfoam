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
    Foam::burstModels::pressure

Description
    burst model using pressure

SourceFiles
    overpressureBurstModel.C
    overpressureBurstModel.C

\*---------------------------------------------------------------------------*/

#ifndef overpressureBurstModel_H
#define overpressureBurstModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "burstModel.H"


namespace Foam
{
namespace burstModels
{

/*---------------------------------------------------------------------------*\
                    Class pressureBurstModel Declaration
\*---------------------------------------------------------------------------*/

class pressure
:
    public burstModel
{
protected:
// Protected variables

    //- Name of pressure field
    word pName_;

    //- Reference pressure
    scalar pRef_;

    //- Pressure causes the patch to burst
    scalar pBurst_;


public:

    //- Runtime type information
    TypeName("pressure");

    // Constructor
    pressure(const dictionary& dict);

    //- Clone
    virtual autoPtr<burstModel> clone() const
    {
        return autoPtr<burstModel>(new pressure(*this));
    }

    //- Destructor
    virtual ~pressure();


    // Member Functions

        //- Return if partial burst is used
        inline const word& pName() const
        {
            return pName_;
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
