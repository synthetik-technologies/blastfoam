/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
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
    Foam::immersedMeshMapper

Description


SourceFiles
    immersedMeshMapper.C

\*---------------------------------------------------------------------------*/

#ifndef immersedMeshMapper_H
#define immersedMeshMapper_H

#include "volFields.H"
#include "RefineBalanceMeshObject.H"
#include "immersedBoundaryObject.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class immersedMeshMapper Declaration
\*---------------------------------------------------------------------------*/

class immersedMeshMapper
:
    public BalanceMeshObject
{
public:
    // Type enumerations

        //- Mesh items to sample
        enum interpolationMethod
        {
            INVERSE_DISTANCE
        };

        static const NamedEnum<interpolationMethod, 1>
            interpolationMethodNames_;

private:
// Private data

    //- Name of object
    const word& name_;

    //- Reference to global polyMesh
    const polyMesh& pMesh_;

    //- Reference to the immersedMesh
    const fvMesh& immersedMesh_;

    //- Reference to the immersed object
    const immersedBoundaryObject* immersedObjectPtr_;

    //- What to sample
    const interpolationMethod mode_;

    word patchName_;

    //- Index of the immersed mesh interface
    label interfaceIndex_;

    //- Interpolation weights from immersed object to immersed mesh
    labelListList immersedObjectToMeshPoints_;

    //- Interpolation weights from immersed object to immersed mesh
    List<scalarList> immersedObjectToMeshWeights_;

    //- Interpolation weights from immersed mesh to immersed object
    labelListList immersedMeshToObjectPoints_;

    //- Interpolation weights from immersed mesh to immersed object
    List<scalarList> immersedMeshToObjectWeights_;

    //- Mapping from immersed mesh to true mesh
    globalIndex* globalImmersedMeshFaceMapPtr_;

    //- Option point displacement field from the immersed mesh
    const pointVectorField* pointDPtr_;
    const volVectorField* DPtr_;


    // Protected Member Functions

        //- Find the nearest face given a set of points
        //  Indexing is based on the interfaceIndex patch
        labelList findNearestFaces(const pointField& points) const;

        //- Calculate mapping
        void calcMapping();


public:

    //- Runtime type information
//     TypeName("immersedMeshMapper");


    // Constructors

    immersedMeshMapper
    (
        const fvMesh& immersedMesh,
        const immersedBoundaryObject& immersedObject,
        const word& patchName
    );

    //- Destructor
    virtual ~immersedMeshMapper();


    // Member Functions

        // Access IBM mapping variable

            //- Return the name of the object
            inline const word& name() const;

            //- Return the immersedMesh
            inline const fvMesh& immersedMesh() const;

            //- Return the immersedObject
            inline const immersedBoundaryObject& immersedObject() const;

            //- Return the index of the interface patch
            inline label interfaceIndex() const;


        // Interpolation and forcing

            //- Interpolate using weights
            template<class Type>
            tmp<Field<Type>> mapImmersedToBoundary(const Field<Type>&) const;

            //- Interpolate using weights
            template<class Type>
            tmp<Field<Type>> mapObjectToBoundary(const Field<Type>&) const;


        //- Set the displacement points
        void setDisplacement(const pointVectorField&, const volVectorField&);

        //- Update mapping
        virtual void updateObject()
        {
            calcMapping();
        }

    // Dummy write
    virtual bool writeData(Ostream&) const
    {
        return true;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "immersedMeshMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "immersedMeshMapperI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
