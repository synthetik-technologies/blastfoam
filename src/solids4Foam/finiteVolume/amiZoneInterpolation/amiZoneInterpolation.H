/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    amiZoneInterpolation

Description
    AMIInterpolation for primitive patches with added point-to-point
    interpolation functions.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    amiZoneInterpolation.C

\*---------------------------------------------------------------------------*/

#ifdef OPENFOAMESIORFOUNDATION

#ifndef amiZoneInterpolation_H
#define amiZoneInterpolation_H

#ifdef OPENFOAMESI
    #include "AMIInterpolation.H"
#else
    #include "newAMIInterpolation.H"
#endif
#include "PrimitivePatch.H"
#include "face.H"
#include "SubList.H"
#include "pointField.H"
#include "PrimitivePatchInterpolation.H"
#include "standAlonePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class amiZoneInterpolation
#ifdef OPENFOAMESI
:
    public AMIInterpolation<standAlonePatch, standAlonePatch>
#else
:
    public newAMIInterpolation<standAlonePatch, standAlonePatch>
#endif
{
    // Private Data

        //- Const reference to source patch
        const standAlonePatch sourcePatch_;

        //- Const reference to target patch
        const standAlonePatch targetPatch_;

        //- Source patch interpolator
        PrimitivePatchInterpolation<standAlonePatch>
            sourcePatchInterp_;

        //- Target patch interpolator
        PrimitivePatchInterpolation<standAlonePatch>
            targetPatchInterp_;

        //- Target face into which each point of source patch is projected
        mutable List<labelPair>* sourcePointAddressingPtr_;

        //- Point weighting factors
        //mutable FieldField<Field, scalar>* sourcePointWeightsPtr_;
        mutable List< List<scalar> >* sourcePointWeightsPtr_;

        //- Distance to intersection for source patch points
        mutable scalarField* sourcePointDistancePtr_;

        //- Source face into which each point of target patch is projected
        mutable List<labelPair>* targetPointAddressingPtr_;

        //- Point weighting factors
        //mutable FieldField<Field, scalar>* targetPointWeightsPtr_;
        mutable List< List<scalar> >* targetPointWeightsPtr_;

        //- Distance to intersection for target patch points
        mutable scalarField* targetPointDistancePtr_;


    // Private Member Functions

        //- Calculate point weights
        void calcSourcePointAddressing() const;

        //- Calculate point weights
        void calcSourcePointWeights() const;

        //- Calculate point weights
        void calcTargetPointAddressing() const;

        //- Calculate point weights
        void calcTargetPointWeights() const;

        //- Disallow default bitwise copy construct
        amiZoneInterpolation(const amiZoneInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const amiZoneInterpolation&);

public:

    //- Runtime type information
    TypeName("amiZoneInterpolation");

    // Constructors

        //- Construct from components
        amiZoneInterpolation
        (
            const standAlonePatch& srcPatch,
            const standAlonePatch& tgtPatch,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool requireMatch = true,
            const interpolationMethod& method = imFaceAreaWeight,
            const scalar lowWeightCorrection = -1,
            const bool reverseTarget = false,
            const bool useGlobalPolyPatch = false
        );

    //- Destructor
    virtual ~amiZoneInterpolation();


    // Member Functions

        // Evaluation

            //- Return reference to point addressing
            const List<labelPair>& sourcePointAddr() const;

            //- Return reference to point weights
            const List< List<scalar> >& sourcePointWeights() const;

            //- Return distance to intersection for patch points
            const scalarField& sourcePointDistanceToIntersection() const;

            //- Return reference to point addressing
            const List<labelPair>& targetPointAddr() const;

            //- Return reference to point weights
            const List< List<scalar> >& targetPointWeights() const;

            //- Return distance to intersection for patch points
            const scalarField& targetPointDistanceToIntersection() const;

            //- Interpolate from source to target point-to-point via faces
            template<class Type>
            tmp<Field<Type>> interpolateToTargetPointsTwoStep
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from target to source point-to-point via faces
            template<class Type>
            tmp<Field<Type>> interpolateToSourcePointsTwoStep
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>::null()
            ) const;

            //- Interpolate from source to target point-to-point directly
            template<class Type>
            tmp<Field<Type> > interpolateToTargetPoints
            (
                const Field<Type>& pf
            ) const;

            //- Interpolate from target to source point-to-point directly
            template<class Type>
            tmp<Field<Type> > interpolateToSourcePoints
            (
                const Field<Type>& pf
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "amiZoneInterpolationTemplates.C"
#endif

#endif

#endif // end of #ifdef OPENFOAMESIORFOUNDATION
// ************************************************************************* //
