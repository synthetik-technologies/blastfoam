/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
08-06-2020  Synthetik Applied Technologies: |   Reload exisiting probe files and
                                                move blastProbes to
                                                the nearest cell/face.
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
    Foam::blastProbes

Description
    Set of locations to sample.

    Call write() to sample and write files.

    If continuing a simulation the old probes files will be trimmed to the
    start time and new values will be appended.


    Example of function object specification:
    \verbatim
    max
    {
        type                blastProbes;
        writeControl        timeStep;
        writeInterval       1;

        probeLocations
        (
            (0 0 0)
            (1 0 0)
            (0 1 0)
            ...
        );
        fields
        (
            p
            T
            rho
        );
        append yes;
        adjustLocations no;
        writeVTK yes;
    }
    \endverbatim

Usage
    \table
        Property          | Description               | Required  | Default
        probeLocations    | List of probe locations   | yes
        fields            | Name of  fields           | yes
        append            | Append to end of old probe files | no | yes
        adjustLocations   | Move blastProbes inside mesh   | no        | no
        writeVTK          | Write the locations a vtk file | no   | no
    \endtable

SourceFiles
    blastProbes.C

\*---------------------------------------------------------------------------*/

#ifndef blastProbes_H
#define blastProbes_H

#include "functionObject.H"
#include "HashPtrTable.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "pointField.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "surfaceMesh.H"
#include "wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class objectRegistry;
class dictionary;
class fvMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                          Class blastProbes Declaration
\*---------------------------------------------------------------------------*/

class blastProbes
:
    public functionObject,
    public pointField
{
protected:

    // Protected classes

        //- Class used for grouping field types
        template<class Type>
        class fieldGroup
        :
            public DynamicList<word>
        {
        public:
            //- Construct null
            fieldGroup()
            :
                DynamicList<word>(0)
            {}
        };

    // Protected member data

        //- Const reference to fvMesh
        const fvMesh& mesh_;

        //- Load fields from files (not from objectRegistry)
        bool loadFromFiles_;

        // Read from dictionary

            //- Names of fields to probe
            wordReList fieldSelection_;

            //- Fixed locations, default = yes
            //  Note: set to false for moving mesh calculations where locations
            //        should move with the mesh
            bool fixedLocations_;

            //- Interpolation scheme name
            //  Note: only possible when fixedLocations_ is true
            word interpolationScheme_;

            //- Are the blastProbes appended to the end of a previous output
            Switch append_;

            //- Switch if update is needed before sampling
            bool needUpdate_;


        // Calculated

            //- Categorized scalar/vector/tensor vol fields
            fieldGroup<scalar> scalarFields_;
            fieldGroup<vector> vectorFields_;
            fieldGroup<sphericalTensor> sphericalTensorFields_;
            fieldGroup<symmTensor> symmTensorFields_;
            fieldGroup<tensor> tensorFields_;

            //- Categorized scalar/vector/tensor surf fields
            fieldGroup<scalar> surfaceScalarFields_;
            fieldGroup<vector> surfaceVectorFields_;
            fieldGroup<sphericalTensor> surfaceSphericalTensorFields_;
            fieldGroup<symmTensor> surfaceSymmTensorFields_;
            fieldGroup<tensor> surfaceTensorFields_;

            //- Cells to be probed (obtained from the locations)
            labelList elementList_;

            //- Physical location of cells to be probed
            //  Used to find nearest point if mesh is refined/unrefined
            pointField elementLocations_;

            //- Faces to be probed
            labelList faceList_;

            //- Current open files
            HashPtrTable<OFstream> probeFilePtrs_;


    // Protected Member Functions

        //- Clear old field groups
        void clearFieldGroups();

        //- Append fieldName to the appropriate group
        label appendFieldGroup(const word& fieldName, const word& fieldType);

        //- Classify field types, returns the number of fields
        label classifyFields();

        //- Return the nearest face to a point
        label findFaceIndex
        (
            const fvMesh& mesh,
            const label celli,
            const vector& pt
        ) const;

        //- Find cells and faces containing blastProbes
        virtual void findElements
        (
            const fvMesh&,
            const bool print,
            const bool movePts = false
        );

        //- Classify field type and Open/close file streams,
        //  returns number of fields to sample
        label prepare();


private:

        //- Sample and write a particular volume field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );


        //- Sample and write a particular surface field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        );

        //- Sample and write all the fields of the given type
        template<class Type>
        void sampleAndWrite(const fieldGroup<Type>&);

        //- Sample and write all the surface fields of the given type
        template<class Type>
        void sampleAndWriteSurfaceFields(const fieldGroup<Type>&);


public:

    //- Runtime type information
    TypeName("blastProbes");


    // Constructors

        //- Construct from Time and dictionary
        blastProbes
        (
            const word& name,
            const Time& time,
            const dictionary& dict
        );

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        blastProbes
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool loadFromFiles = false
        );

        //- Disallow default bitwise copy construction
        blastProbes(const blastProbes&) = delete;


    //- Destructor
    virtual ~blastProbes();


    // Member Functions

        //- Return names of fields to probe
        virtual const wordReList& fieldNames() const
        {
            return fieldSelection_;
        }

        //- Return locations to probe
        virtual const pointField& probeLocations() const
        {
            return *this;
        }

        //- Return location for probe i
        virtual const point& probe(const label i) const
        {
            return operator[](i);
        }

        //- Cells to be probed (obtained from the locations)
        const labelList& elements() const
        {
            return elementList_;
        }

        //- Read the blastProbes
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&);

        //- Update for changes of mesh due to readUpdate
        virtual void readUpdate(const polyMesh::readUpdateState state)
        {}

        //- Sample a volume field at all locations
        template<class Type>
        tmp<Field<Type>> sample
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Sample a single vol field on all sample locations
        template<class Type>
        tmp<Field<Type>> sample(const word& fieldName) const;

        //- Sample a single scalar field on all sample locations
        template<class Type>
        tmp<Field<Type>> sampleSurfaceFields(const word& fieldName) const;

        //- Sample a surface field at all locations
        template<class Type>
        tmp<Field<Type>> sample
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const blastProbes&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "blastProbesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
