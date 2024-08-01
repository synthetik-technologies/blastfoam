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
    Foam::displacementConstraint

Description
    Base class for constraining nodal displacement and velocity

SourceFiles
    displacementConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef displacementConstraint_H
#define displacementConstraint_H

#include "pointFields.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class displacementConstraint
{
public:
    // Read a list of component types i.e. (x y) or (xx xy yz);
    static List<label> readComponents(Istream& is);

    enum selectionType
    {
        UNKNOWN = -1,
        ZONE,
        SET,
        PATCH
    };

    static const Foam::NamedEnum<selectionType, 3> selectionTypeNames;

protected:

    //- Name of the constraint
    word name_;

    //- Reference to the polyMesh
    const polyMesh& mesh_;

    //- Displacement field
    pointVectorField& D_;

    //- Velocity field
    pointVectorField& U_;

    //- Selection type
    selectionType selectionType_;

    //- Name of the selection
    word selectionName_;

    //- Constrained directions
    List<label> dims_;

    //- Constrained nodes
    List<label> nodes_;


    //- Update constrained nodes
    void updateNodes();


public:

    //- Runtime type information
    TypeName("displacementConstraint")

        // Declare runtime constructor selection tables
        declareRunTimeSelectionTable
        (
            autoPtr,
            displacementConstraint,
            dictionary,
            (
                const word& name,
                pointVectorField& D,
                pointVectorField& U,
                const dictionary& dict
            ),
            (name, D, U, dict)
        );

        // Constructor
        displacementConstraint
        (
            const word& name,
            pointVectorField& D,
            pointVectorField& U,
            const dictionary& dict
        );

        //- Destructor
        virtual ~displacementConstraint();

        // Selector
        static autoPtr<displacementConstraint> New
        (
            pointVectorField& D,
            pointVectorField& U,
            const dictionary& dict
        );

        // Selector
        static autoPtr<displacementConstraint> New
        (
            const word& type,
            const word& name,
            pointVectorField& D,
            pointVectorField& U,
            const dictionary& dict
        );

        //- Class used for the read-construction of
        //  PtrLists of displacementConstraint.
        class iNew
        {
            pointVectorField& D_;
            pointVectorField& U_;

        public:

            iNew
            (
                pointVectorField& D,
                pointVectorField& U
            )
            :
                D_(D),
                U_(U)
            {}

            autoPtr<displacementConstraint> operator()(Istream& is) const
            {
//                 word constraintType(is);
                word constraintName(is);
                dictionary dict(is);
                word constraintType(dict.lookup("type"));

                return displacementConstraint::New
                (
                    constraintType,
                    constraintName,
                    D_,
                    U_,
                    dict
                );
            }
        };

        //- Clone
        virtual autoPtr<displacementConstraint> clone() const = 0;


    // Memeber Functions

        // Access

            //- Return the name of the constraint
            const word& name() const
            {
                return name_;
            }

            //- Return the constrained nodes
            const labelList& nodes() const
            {
                return nodes_;
            }

            //- Evaluate the coefficient at the given point
            virtual void constrain() = 0;
};



} // End namespace Foam

#endif
