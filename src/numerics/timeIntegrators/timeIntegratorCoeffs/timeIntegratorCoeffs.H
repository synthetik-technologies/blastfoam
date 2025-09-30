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
    Foam::timeIntegratorCoeffs

Description
    Base class for time integration coefficients

SourceFiles
    timeIntegratorCoeffs.C
    timeIntegratorCoeffsNew.C

\*---------------------------------------------------------------------------*/

#ifndef timeIntegratorCoeffs_H
#define timeIntegratorCoeffs_H

#include "runTimeSelectionTables.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class timeIntegratorCoeffs Declaration
\*---------------------------------------------------------------------------*/

class timeIntegratorCoeffs
{
protected:
    // Protected data

        //- Number of steps
        label nSteps_;


public:

    //- Runtime type information
    TypeName("timeIntegratorCoeffs");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            timeIntegratorCoeffs,
            dictionary,
            (Istream& is),
            (is)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            timeIntegratorCoeffs,
            dictionaryEmbedded,
            (Istream& is),
            (is)
        );

    // Constructor
    timeIntegratorCoeffs(const label nSteps);


    //- Destructor
    virtual ~timeIntegratorCoeffs();


    // Selectors

        //- Return RK coefficnents
        static autoPtr<timeIntegratorCoeffs> New
        (
            const dictionary& dict,
            const bool embedded
        );


    // Member Functions

        // Coefficient functions

            //- Return the number of steps
            inline label nSteps() const
            {
                return nSteps_;
            }

            //- Set coefficients
            virtual void set
            (
                List<List<scalar>>& as,
                List<List<scalar>>& bs,
                const label index
            ) const = 0;

            //- Update coefficients
            virtual bool update(const label index, const scalar t);

            //- Return the index of saved old fields
            virtual List<label> oldIs(const List<List<scalar>>& as) const;

            //- Return the index of saved delta fields
            virtual List<label> deltaIs(const List<List<scalar>>& bs) const;

            // Set the initial and end time coefficients (all based on current time)
            virtual void setTimeFactors
            (
                const List<List<scalar>>& as,
                const List<List<scalar>>& bs,
                List<scalar>& f0,
                List<scalar>& f
            ) const;

            //- Return max Courant number
            virtual scalar maxCo() const
            {
                return 1.0;
            }



        // Mapping/save functions

            //- Are there any old times that need to be saved
            virtual bool save() const
            {
                return false;
            }

            //- Return the map for saving old values
            virtual List<label> oldSaveMap() const
            {
                return List<label>();
            }

            //- Return the map for saving old values
            virtual List<label> deltaSaveMap() const
            {
                return List<label>();
            }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
