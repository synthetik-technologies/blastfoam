/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::boundaryTraction

Description
    Base class for appling boundary tractions/loads

SourceFiles
    boundaryTraction.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryTraction_H
#define boundaryTraction_H

#include "feMesh1.H"
#include "polyMesh.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class boundaryTraction
{
protected:

    //- Name of the boundary traction
    const word name_;

    //- Reference to the femesh1
    const feMesh1& femesh_;

    //- Reference to polyMesh
    const polyMesh& mesh_;

    //- Optional pointer to displacement field
    const pointVectorField* DPtr_;

    //- Patch index
    label patchID_;

    //- Update elements
    void updateElements();

public:

    //- Runtime type information
    TypeName("boundaryTraction")

        // Declare runtime constructor selection tables
        declareRunTimeSelectionTable
        (
            autoPtr,
            boundaryTraction,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const feMesh1& femesh,
                const pointVectorField* DPtr
            ),
            (name, dict, femesh, DPtr)
        );

        // Constructor
        boundaryTraction
        (
            const word& name,
            const dictionary& dict,
            const feMesh1& femesh,
            const pointVectorField* DPtr
        );

        //- Destructor
        virtual ~boundaryTraction();

        // Selector
        static autoPtr<boundaryTraction> New
        (
            const word& name,
            const dictionary& dict,
            const feMesh1& femesh,
            const pointVectorField* DPtr
        );

        //- Class used for the read-construction of
        //  PtrLists of boundaryTraction.
        class iNew
        {
            const feMesh1& femesh_;
            const pointVectorField* DPtr_;

        public:

            iNew(const feMesh1& femesh, const pointVectorField* DPtr)
            :
                femesh_(femesh),
                DPtr_(DPtr)
            {}

            autoPtr<boundaryTraction> operator()(Istream& is) const
            {
                const word name(is);
                dictionary dict(is);

                return boundaryTraction::New
                (
                    name,
                    dict,
                    femesh_,
                    DPtr_
                );
            }
        };

        //- Clone
        virtual autoPtr<boundaryTraction> clone() const = 0;


    // Member Functions

        // Access

            //- Return patchIndex
            const label& patchID() const
            {
                return patchID_;
            }


        // Edit

            //- Update values
            virtual void update() = 0;

            //- Evaluate the coefficient at the given point
            virtual vector getTraction
            (
                const labelList& nodeLabels,
                const scalarList& shape,
                const vector& n
            ) const;

            // Return pressure given shape function and nodes
            virtual scalar pressure
            (
                const labelList& nodeLabels,
                const scalarList& shape,
                const vector& n
            ) const = 0;

            // Return traction given shape function, nodes, and normal vector
            virtual vector traction
            (
                const labelList& nodeLabels,
                const scalarList& shape,
                const vector& n
            ) const = 0;
};



} // End namespace Foam

#endif
