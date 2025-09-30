/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    externalPressureCellRemovalLaw

Description
    Cell removal law which selects boundary cells that meet some external
    condition to remove.

    The removed cells are transferred to a subMesh and the map to the external
    mesh is calculated using the meshToMesh class

SourceFiles
    externalConditionCellRemovalLaw.C

\*---------------------------------------------------------------------------*/

#ifndef externalPressureCellRemovalLaw_H
#define externalPressureCellRemovalLaw_H

#include "cellRemovalLaw.H"
#include "fvMeshSubset.H"
#include "meshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class externalConditionCellRemovalLaw Declaration
\*---------------------------------------------------------------------------*/

class externalConditionCellRemovalLaw
:
    public cellRemovalLaw
{
    // Private Data Member

        // Critical value of pressure
        const scalar pCritical_;

        // Name of the pressure field
        const word pName_;

        //- Name of coupled patch
        const label coupledPatchID_;

        // Index of patch to put exposed faces in
        const label patchID_;

        //- Removed subMesh
        mutable autoPtr<fvMeshSubset> subsetter_;

        //- MeshToMesh mapping
        mutable autoPtr<meshToMesh> meshMap_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        externalConditionCellRemovalLaw(const externalConditionCellRemovalLaw&);

        //- Disallow default bitwise assignment
        void operator=(const externalConditionCellRemovalLaw&);

        //- Save geometric fields
        template<class Type, template<class> class Patch, class Mesh>
        void saveGeoFields(objectRegistry&) const;


public:

    //- Runtime type information
    TypeName("externalCondition");

    // Static data members


    // Constructors

        //- Construct from dictionary
        externalConditionCellRemovalLaw
        (
            const word& name,
            fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~externalConditionCellRemovalLaw();


    // Member Functions

        // Return list of cell IDs to be removed
        virtual labelList cellsToRemove();

        // Return patch index for the faces to be exposed
        virtual label exposedFacesPatchID();

        //- Return if cells have been removed
        bool active() const;

        //- Return the meshToMesh mapping
        const meshToMesh& meshMap() const;

        //- Interpolate field and combine
        template<class Type, class CombineOp>
        void map
        (
            const word&,
            const CombineOp&,
            DimensionedField<Type, volMesh>&
        ) const;
};


template<>
void externalConditionCellRemovalLaw::map
(
    const word&,
    const plusEqOp<scalar>&,
    DimensionedField<scalar, volMesh>&
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "externalConditionCellRemovalLawTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
