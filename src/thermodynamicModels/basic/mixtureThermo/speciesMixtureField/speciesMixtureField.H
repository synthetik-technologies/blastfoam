/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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
    Template class to define a runtime selectable fluidThermo model

SourceFiles
    speciesMixtureField.C

\*---------------------------------------------------------------------------*/

#ifndef speciesMixtureField_H
#define speciesMixtureField_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicThermoModel.H"
#include "volFields.H"
#include "MeshObject.H"
#include "PtrList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class speciesMixtureField Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class speciesMixtureField
:
    public PtrList<ThermoType>,
    public UpdateableMeshObject<polyMesh>
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
    PtrList<PtrList<ThermoType>> faceMixtures_;


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

        //- Update cell mixtures
        void updateCellMixtures();

        //- Update face mixtures for a given patch
        void updateFaceMixtures(const label patchi);

        //- Update all mixtures
        void update();

        //- Update size of fields
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update size of fields
        virtual bool movePoints()
        {
            return true;
        }


    // Access functions

        //- Return the boundary mixtures
        const PtrList<PtrList<ThermoType>>& boundary() const
        {
            return faceMixtures_;
        }

        //- Return the boundary mixtures
        PtrList<PtrList<ThermoType>>& boundary()
        {
            return faceMixtures_;
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
