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
#include "patchToPatch.H"
#include "fvPatch.H"


namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class burstModel Declaration
\*---------------------------------------------------------------------------*/

class burstModel
{
protected:
// Protected variables

    //- Constructing dictionary
    const dictionary dict_;

    //- Is this for a coupled set of patches
    word patchToPatchType_;

    //- Is a partial burst used
    bool partialBurst_;

    //- Is the delta across coupled patches used
    bool useDelta_;

    //- Use average pressure to determine burst (only for non partial burst)
    bool useAverage_;

    //- Log information
    bool log_;

    //- Is an update needed
    mutable bool needUpdate_;

    //- Patch to patch storage
    mutable autoPtr<patchToPatch> mappingPtr_;

    // Protected member function

        //- Helper function to find faces that have burst
        bool findBurstFaces
        (
            const scalarField& pf1,
            const fvPatch&  patch1,
            const scalarField& pf2,
            const fvPatch&  patch2,
            const scalar burstValue,
            labelHashSet& master,
            labelHashSet& slave
        ) const;

public:

    //- Runtime type information
    TypeName("burstModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            burstModel,
            dictionary,
            (
                const dictionary& dict,
                const bool coupled
            ),
            (dict, coupled)
        );

    //- Constructor
    burstModel
    (
        const dictionary& dict,
        const bool coupled
    );


    //- Destructor
    virtual ~burstModel();


    // Selectors

        //- Return a new burstModel object
        static autoPtr<burstModel> New
        (
            const dictionary& dict,
            const bool coupled
        );

        //- Clone
        virtual autoPtr<burstModel> clone() const = 0;


    // Member Functions

        //- Set if an update is needed
        void needUpdate() const;

        //- Update the patch to patch mapping
        void updateMapping
        (
            const primitivePatch& patch1,
            const primitivePatch& patch2
        ) const;

        //- Find and set face indices to change, return if any added
        virtual bool facesToChange
        (
            const fvPatch& patch1,
            const fvPatch& patch2,
            labelList& masterFaces,
            labelList& slaveFaces
        ) const = 0;

        //- Find and set face indices to change, return if any added
        virtual bool facesToChange
        (
            const fvPatch& patch,
            labelList& faces
        ) const = 0;

        //- Return if burst occurs
        //  1 owner side
        //  -1 neighbour side
        // 0 no burst
        virtual label update
        (
            const objectRegistry& obr,
            const labelList& ownCells,
            const labelList& neiCells,
            const scalarField& V
        ) const = 0;

        //- Write
        virtual void writeData(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
