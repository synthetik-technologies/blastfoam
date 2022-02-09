/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-12-11 Jeff Heylmun:    Added storage of mixture fields and derivation
                            from UpdateableMeshObject
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
    Foam::speciesMixtureField

Description
    Basic species mixture class responsible for storing and updating mixtures.
    Derived from UpdateableMeshObject so checking of sizes does not need to be
    done, and the updating of sizes is done in updateMesh.

SourceFiles
    speciesMixtureField.C

\*---------------------------------------------------------------------------*/

#ifndef speciesMixtureField_H
#define speciesMixtureField_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "RefineBalanceMeshObject.H"
#include "PtrList.H"
#include "emptyPolyPatch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class speciesMixtureField Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class speciesMixtureField
:
    public PtrList<ThermoType>,
    public MeshObject
    <
        fvMesh,
        DistributeableMeshObject,
        speciesMixtureField<ThermoType>
    >

{
// Private data

    //- Const reference to mesh
    const fvMesh& mesh_;

    //- Const reference to mass fractions
    const PtrList<volScalarField>& Ys_;

    //- Species data
    const PtrList<ThermoType>& speciesData_;

    //- Temporary storage for the cell/face mixture thermo data
    mutable ThermoType mixture_;

    //- Stored thermo state at boundaries
    PtrList<ThermoType> faceMixtures_;

    //- Start indexes of patches
    labelList start_;

    //- Return label index
    label patchFaceIndex(const label patchi, const label facei) const
    {
        return start_[patchi] + facei;
    }

    //- Return the number of boundary faces
    label nBoundaryFaces() const
    {
        label nf = 0;
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        forAll(pbm, patchi)
        {
            if (!isA<emptyPolyPatch>(pbm[patchi]))
            {
                nf += pbm[patchi].size();
            }
        }
        return nf;
    }

    //- Update size of fields
    void updateCells();
    void updatePatches(const bool validBoundary);


public:

    // Constructor
    speciesMixtureField
    (
        const fvMesh& mesh,
        const PtrList<volScalarField>& Ys,
        const PtrList<ThermoType>& speciesData,
        const word& phaseName
    );

    //- Destructor
    virtual ~speciesMixtureField();


    // Member Functions

        //- Update cell mixture for celli
        const ThermoType& cellMixture(const label celli) const;

        //- Update mixture for a given face on a patch
        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        //- Update all mixtures
        void updateMixture();


    // Access functions

        //- Return the boundary mixtures
        const ThermoType& boundary
        (
            const label patchi,
            const label facei
        ) const
        {
            return faceMixtures_[patchFaceIndex(patchi, facei)];
        }

        //- Return the boundary mixtures
        ThermoType& boundary
        (
            const label patchi,
            const label facei
        )
        {
            return faceMixtures_[patchFaceIndex(patchi, facei)];
        }

        //- List operators
        const ThermoType& operator[](const label i) const
        {
            return PtrList<ThermoType>::operator[](i);
        }

        ThermoType& operator[](const label i)
        {
            return PtrList<ThermoType>::operator[](i);
        }

        virtual bool writeData(Ostream&) const
        {
            return true;
        }

    // Mesh object functions

        virtual bool movePoints()
        {
            return false;
        }

        virtual void reorderPatches
        (
            const labelUList& newToOld,
            const bool validBoundary
        )
        {
            updatePatches(validBoundary);
        }
        virtual void addPatch(const label patchi)
        {
            updatePatches(true);
        }

        //- Save the intact field so it can be mapped after the mesh is distributed
        virtual void preDistribute()
        {}

        //- Distribute the intact field
        virtual void distribute(const mapDistributePolyMesh&)
        {
            updateCells();
            updatePatches(true);
        }

        //- Update the intact field after mapping
        virtual void updateMesh(const mapPolyMesh& mpm)
        {
            updateCells();
            updatePatches(true);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "speciesMixtureField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
