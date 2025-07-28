/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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
    Foam::MeshField

SourceFiles
    MeshField.C
    MeshFields.C

Description
    Regiseterd field that is updated by the mesh

\*---------------------------------------------------------------------------*/

#ifndef MeshField_H
#define MeshField_H

#include "IOField.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class MeshField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class MeshField
:
    public IOField<Type>
{
public:
    // Public typedefs

    typedef Field<Type> FieldType;
    typedef typename GeoMesh::Mesh Mesh;

private:
    // Private member functions

        //- Reference to the mesh
        const polyMesh& mesh_;

        //- Ensure a given field is the correct size
        void checkSize(const Field<Type>& f) const
        {
            if (f.size() != GeoMesh::size(this->mesh_))
            {
                FatalErrorInFunction
                    << "size of field = " << f.size()
                    << " is not the same as the size of mesh = "
                    << GeoMesh::size(this->mesh_)
                    << abort(FatalError);
            }
        }

public:

    TypeName("MeshField");

    // Constructors

        //- Construct from IOobject
        MeshField(const IOobject&, const polyMesh&);

        //- Construct from IOobject and size (does not set values)
        MeshField(const IOobject&, const polyMesh&, const Type& val);

        //- Construct from components
        MeshField(const IOobject&, const polyMesh&, const Field<Type>&);

        //- Move construct by transferring the Field contents
        MeshField(const IOobject&, const polyMesh&, Field<Type>&&);

        //- Move construct by transferring the Field contents
        MeshField(const IOobject&, const polyMesh&, const tmp<Field<Type>>&);


    //- Destructor
    virtual ~MeshField();


    // Member functions

        //- Return the mesh
        inline const polyMesh& mesh() const
        {
            return mesh_;
        }

    // Member operators
    using Field<Type>::operator=;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MeshField.C"
#endif

#endif

// ************************************************************************* //
