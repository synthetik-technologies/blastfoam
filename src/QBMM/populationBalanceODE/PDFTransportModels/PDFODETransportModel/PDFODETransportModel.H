/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 Alberto Passalacqua
     \\/     M anipulation  |
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
    Foam::PDFODETransportModel

Description
    Abstract class for the run-time selection of the PDF transport model.

SourceFiles
    PDFODETransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef PDFODETransportModel_H
#define PDFODETransportModel_H

#include "PDFTransportModel.H"
#include "fvMesh.H"
#include "volFields.H"
#include "integrationSystem.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PDFODETransportModel Declaration
\*---------------------------------------------------------------------------*/

class PDFODETransportModel
:
    public integrationSystem
{
    // Private data

         //- Name of the PDFODETransportModel
        const word name_;


protected:

    // Protected data

        //- Constant reference to the mesh
        const fvMesh& mesh_;


    // Protected member functions

        //- Return number of moments
        virtual label nMoments() const = 0;


public:

    // Constructors

        //- Construct from components
        PDFODETransportModel
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construct
        PDFODETransportModel(const PDFODETransportModel&) = delete;


    //- Destructor
    virtual ~PDFODETransportModel();


    // Member Functions


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const PDFODETransportModel&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
