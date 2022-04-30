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
    Foam::burstModel

Description
    Base class for time integration

SourceFiles
    burstModel.C
    burstModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef burstModel_H
#define burstModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "runTimeSelectionTables.H"
#include "fvPatchFields.H"


namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class burstModel Declaration
\*---------------------------------------------------------------------------*/

class burstModel
{
protected:
// Protected variables

    //- Is a partial burst used
    bool partialBurst_;

    //- Has the patch completely burst
    mutable bool burst_;


public:

    //- Runtime type information
    TypeName("burstModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            burstModel,
            dictionary,
            (const dictionary& dict),
            (dict)
        );

    // Constructor
    burstModel(const dictionary& dict);
    burstModel();


    //- Destructor
    virtual ~burstModel();


    // Selectors

        //- Return a new burstModel object
        static autoPtr<burstModel> New(const dictionary& dict);

        //- Clone
        virtual autoPtr<burstModel> clone() const = 0;


    // Member Functions

        //- Return if partial burst is used
        inline bool partialBurst() const
        {
            return partialBurst_;
        }


        //- Update the given intact field
        virtual bool update
        (
            const fvPatch&,
            scalarField& intact
        ) const = 0;

        //- Write
        virtual void writeData(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
