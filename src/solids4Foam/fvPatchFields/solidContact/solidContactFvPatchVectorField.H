/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    solidContactFvPatchVectorField

Description
    Contact boundary conditions for finite volume solidMechanics solvers.

    The treatment of the normal contact and tangential contact can be specified
    individually at run-time using the normalContactModel and
    frictionContactModel.

    The boundary condition keeps a copy of the master and slave patch and
    keeps them in the deformed configuration; this allows the calculation of
    contact distances and interpolation between the globalPatchZones.

    The distance calculations and interpolations are performed by the GGI class.

    More details in:

    P. Cardiff, A. Karać, A. Ivanković: Development of a Finite Volume contact
    solver based on the penalty method. Computational Materials Science,
    03/2014, 64:283–284. DOI: 10.1016/j.commatsci.2012.03.011.

    P. Cardiff, Z. Tukovic, P. De Jaeger, M. Clancy and A. Ivankovic. A
    Lagrangian cell-centred finite volume method for metal forming simulation,
    doi=10.1002/nme.5345.

    P. Cardiff, Development of the Finite Volume Method for Hip Joint Stress
    Analysis, PhD thesis, University College Dublin, 2012.


SourceFiles
    solidContactFvPatchVectorField.C
    solidContactFvPatchVectorFieldCalcContact.C

Author
    Philip Cardiff, UCD. All rights reserved.
    Rewrite by Hrvoje Jasak.

\*---------------------------------------------------------------------------*/

#ifndef solidContactFvPatchVectorField_H
#define solidContactFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidTractionFvPatchVectorField.H"
#include "normalContactModel.H"
#include "frictionContactModel.H"
#include "globalPolyPatch.H"
#include "ggiInterpolation.H"
#include "Switch.H"
#include "boundBox.H"
#include "standAlonePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class solidContactFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidContactFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{
    // Private data

        //- Store a copy of the patch dictionary
        //  This allows us to use lazy evaluation when creating the contact
        //  models
        const dictionary dict_;

        //- Is this patch the master
        const Switch master_;

        //- Write the deformed contact zones to VTK surface files
        const Switch writeZoneVTK_;

        //- Write the contact point distances fields for viewing
        const Switch writePointDistanceFields_;

        //- Shadow patch names
        mutable wordList shadowPatchNames_;

        //- Shadow patch indices
        mutable labelList* shadowPatchIndicesPtr_;

        //- Is the contact boundary active
        Switch rigidMaster_;

        //- Normal contact model pointers
        mutable PtrList<normalContactModel> normalModels_;

        //- Friction contact model pointer
        mutable PtrList<frictionContactModel> frictionModels_;

        //- Store current penalty scale factors (if using penalty contact model)
        mutable scalarField normalPenaltyFactors_;

        //- Master face zone kept in deformed configuration
        mutable globalPolyPatch* zonePtr_;

        //- Shadow face zones kept in deformed configuration
        //  Stored by the master patch
        mutable PtrList<globalPolyPatch> shadowZones_;

        //- Zone-to-zone interpolations
        mutable PtrList<ggiStandAlonePatchInterpolation> zoneToZones_;

        //- GGI quick reject algorithm
        ggiInterpolation::quickReject quickReject_;

        //- Top corner of the region where te contact is expected to be
        // This aims to save time in the GGI distance calculation
        const vector regionOfInterestTopCorner_;

        //- Bottom corner of the region where te contact is expected to be
        // This aims to save time in the GGI distance calculation
        const vector regionOfInterestBottomCorner_;

        //- Region where the contact is expected to be
        //  This aims to save time in the GGI distance calculations
        const boundBox regionOfInterest_;

        //- Contact field for all the contacts, stored on the current patch
        //  1 means in contact
        //  0 means not in contact
        scalarField contact_;

        //- Contact field corresponding to each shadow patch, stored on the
        //  current patch
        //  1 means in contact
        //  0 means not in contact
        mutable PtrList<scalarField> contactPerShadow_;

        //- Flag to enable approach for scaling traction near the downstream
        //  patch to fix an issue with unphysical behaviour
        Switch scaleFaceTractionsNearDownstreamPatch_;

        //- Field to scale traction of slave patch
        mutable autoPtr<scalarField> scaleTractionFieldPtr_;

        //- Current time index
        label curTimeIndex_;


    // Private Member Functions

        //- Is a moving mesh (e.g. updated Lagrangian) approach be used
        bool movingMesh() const;

        //- Make the shadow patch names
        void makeShadowPatchNames(const dictionary& dict) const;

        //- Calculate the shadow patch indices
        void calcShadowPatchIndices() const;

        //- Make normal contact models
        void makeNormalModels(const dictionary& dict) const;

        //- Make friction contact models
        void makeFrictionModels(const dictionary& dict) const;

        //- Clear out demand driven data
        //virtual void clearOut();

        //- Create the zone
        void calcZone() const;

        //- Create the shadow zones
        void calcShadowZones() const;

        //- Create the zone-to-zone interpolators
        void calcZoneToZones() const;

        //- Move the master and slave zones to the deformed configuration
        void moveZonesToDeformedConfiguration();

        // Set the contactPerShadow field
        void calcContactPerShadow() const;

        // //- Interpolate zone face values to points
        // template<class Type>
        // tmp<Field<Type> > zoneFaceToPointInterpolate
        // (
        //     const label zoneIndex,
        //     const Field<Type>& zoneField
        // ) const;

        // //- Interpolate zone point values to faces
        // template<class Type>
        // tmp<Field<Type> > zonePointToFaceInterpolate
        // (
        //     const label zoneIndex,
        //     const Field<Type>& zonePointField
        // ) const;

        //- Return field to scale the traction field on the slave patch
        scalarField scaleTractionField() const;

        //- Make the scaleTractionField
        void makeScaleTractionField() const;


public:

    //- Runtime type information
    TypeName("solidContact");


    // Constructors

        //- Construct from patch and internal field
        solidContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidContactFvPatchVectorField onto
        //  a new patch
        solidContactFvPatchVectorField
        (
            const solidContactFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidContactFvPatchVectorField
        (
            const solidContactFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        solidContactFvPatchVectorField
        (
            const solidContactFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new solidContactFvPatchVectorField(*this, iF)
            );
        }

    //- Destructor
    virtual ~solidContactFvPatchVectorField();


    // Member functions

        // Access

            //- Return if this patch is the slave
            Switch master() const
            {
                return master_;
            }

            //- Return the shadow patch names
            const wordList& shadowPatchNames() const;

            //- Return the shadow patch indices
            const labelList& shadowPatchIndices() const;

            //- Return the shadow patch
            //  This function is not allowed if shadowPatchIndices.size() > 1
            const solidContactFvPatchVectorField& shadowPatchField() const;

            //- Return reference to the normal contact models
            PtrList<normalContactModel>& normalModels();

            //- Return const reference to the normal contact models
            const PtrList<normalContactModel>& normalModels() const;

            //- Return reference to the friction contact models
            PtrList<frictionContactModel>& frictionModels();

            //- Return const reference to the friction contact models
            const PtrList<frictionContactModel>& frictionModels() const;

            //- Return normal contact model for the slave patch calling this
            //  function
            normalContactModel& normalModelForThisSlave();

            //- Return friction contact model for the slave patch calling this
            //  function
            frictionContactModel& frictionModelForThisSlave();

            //- Return const reference to deformed zone
            const globalPolyPatch& zone() const;

            //- Return reference to deformed zone
            globalPolyPatch& zone();

            //- Return const reference to deformed shadow zones
            const PtrList<globalPolyPatch>& shadowZones() const;

            //- Return reference to deformed shadow zones
            PtrList<globalPolyPatch>& shadowZones();

            //- Return const reference to patch-to-patch interpolation
            const PtrList<ggiStandAlonePatchInterpolation>&
                zoneToZones() const;

            //- Return reference to patch-to-patch interpolation
            PtrList<ggiStandAlonePatchInterpolation>& zoneToZones();

            //- Return const reference to zone-to-zone interpolation for the
            //  current slave patch; an error is thrown if the master calls this
            //  function
            const ggiStandAlonePatchInterpolation&
            zoneToZoneForThisSlave() const;

            //- Return the global patch for the current slave patch
            //  An error is thrown if the master calls this function
            const globalPolyPatch& zoneForThisSlave() const;

            //- Return reference to contact field
            const scalarField& contact() const
            {
                return contact_;
            }

            //- Return reference to contactPerShadow fields
            PtrList<scalarField>& contactPerShadow();

            //- Return const reference to contactPerShadow fields
            const PtrList<scalarField>& contactPerShadow() const;


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // Member functions

            //- Clear out demand driven data
            virtual void clearOut();

            //- Return the rate of dissipated energy due to friction for the
            //  master patch
            virtual tmp<scalarField> frictionHeatRate() const;

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

//#ifdef NoRepository
//#   include "solidContactFvPatchVectorFieldTemplates.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
