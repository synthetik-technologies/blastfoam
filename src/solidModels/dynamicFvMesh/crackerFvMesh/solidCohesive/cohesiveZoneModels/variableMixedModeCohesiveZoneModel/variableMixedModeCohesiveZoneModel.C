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

\*---------------------------------------------------------------------------*/

#include "variableMixedModeCohesiveZoneModel.H"
#include "addToRunTimeSelectionTable.H"
#include "solidCohesiveFvPatchVectorField.H"
#include "cohesiveZoneInitiation.H"
#include "directFvPatchFieldMapper.H"
#include "crackerFvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(variableMixedModeCohesiveZoneModel, 0)
    addToRunTimeSelectionTable
    (
        cohesiveZoneModel, variableMixedModeCohesiveZoneModel, dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::vector Foam::variableMixedModeCohesiveZoneModel::damageTractionN
(
    const scalar faceDeltaN,
    const scalar faceDeltaS,
    const scalar faceSigmaMax,
    const scalar faceTauMax,
    const vector& faceN
) const
{
    return
        faceN*faceSigmaMax*faceDeltaN
       /(
            Foam::sqrt
            (
                pow(faceDeltaN, 2)
              + pow(faceDeltaS, 2)*pow(faceSigmaMax, 2)/pow(faceTauMax, 2)
            )
        );
}


Foam::vector Foam::variableMixedModeCohesiveZoneModel::damageTractionS
(
    const scalar faceDeltaN,
    const scalar faceDeltaS,
    const scalar faceSigmaMax,
    const scalar faceTauMax,
    const vector& faceN,
    const vector& faceDelta
) const
{
    if (mag(faceDeltaS) > SMALL)
    {
        scalar magTracS =
            faceTauMax*faceDeltaS
            /(
                Foam::sqrt
                (
                    pow(faceDeltaS, 2)
                    + pow(faceDeltaN, 2)
                    *pow(faceTauMax, 2)/pow(faceSigmaMax, 2)
                )
            );

        // Shear faceDelta direction
        const vector slip = ((I - sqr(faceN)) & faceDelta);

        // Slip magnitude
        scalar magSlip = mag(slip);

        // Slip direction
        const vector sDir = slip/magSlip;

        // Return shear traction vector
        return magTracS*sDir;
    }

    // No slip: return zero traction
    return vector::zero;
}


void Foam::variableMixedModeCohesiveZoneModel::calcPenaltyFactor() const
{
    if (patch().size())
    {
        // Calculate penalty factor similar to standardPenalty contact model
        // approx penaltyFactor from mechanical properties
        // this can then be scaled using the penaltyScale

        const label patchID = patch().index();
        const fvMesh& mesh = this->mesh();

        // Lookup the implicit stiffness as a measure of the mechanical
        // stiffness
        const scalar impK =
            gAverage
            (
                mesh.lookupObject<volScalarField>
                (
                    "impK"
                ).boundaryField()[patchID]
            );

        // Average contact patch face area
        scalar faceArea = gAverage(patch().magSf());

        // average contact patch cell volume
        scalar cellVolume = 0.0;

        const volScalarField::Internal& V = mesh.V();
        const unallocLabelList& faceCells =
            mesh.boundary()[patchID].faceCells();

        forAll(mesh.boundary()[patchID], facei)
        {
            cellVolume += V[faceCells[facei]];
        }
        cellVolume /= patch().size();

        // Approximate penalty factor based on:
        // Hallquist, Goudreau, Benson - 1985 - Sliding interfaces with
        // contact-impact in large-scale Lagrangian computations, Computer
        // Methods in Applied Mechanics and Engineering, 51:107-137.
        // We approximate penalty factor for traction instead of force, so we
        // divide their proposed penalty factor by face-area
        penaltyFactor_ = penaltyScale_*impK*faceArea/cellVolume;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from dictionary
Foam::variableMixedModeCohesiveZoneModel::variableMixedModeCohesiveZoneModel
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict
)
:
    cohesiveZoneModel(name, patch, dict),
    sigmaMax_("sigmaMax", dimPressure, dict),
    tauMax_("tauMax", dimPressure, dict),
    GIc_("GIc", dimensionSet(1, 0, -2, 0, 0, 0, 0), dict),
    GIIc_("GIIc", dimensionSet(1, 0, -2, 0, 0, 0, 0), dict),
    cracked_(patch.size(), false),
    tractionN_(patch.size(), 0.0),
    oldTractionN_(patch.size(), 0.0),
    tractionS_(patch.size(), 0.0),
    oldTractionS_(patch.size(), 0.0),
    deltaN_(patch.size(), 0.0),
    oldDeltaN_(patch.size(), 0.0),
    deltaS_(patch.size(), 0.0),
    oldDeltaS_(patch.size(), 0.0),
    deltaEff_(patch.size(), 0.0),
    unloadingDeltaEff_(patch.size(), 0.0),
    GI_(patch.size(), 0.0),
    oldGI_(patch.size(), 0.0),
    GII_(patch.size(), 0.0),
    oldGII_(patch.size(), 0.0),
    penaltyScale_(dict.lookupOrDefault<scalar>("penaltyScale", 1.0)),
    penaltyFactor_(dict.lookupOrDefault<scalar>("penaltyFactor", 0.0))
{
    if (dict.found("cracked"))
    {
        cracked_ = Field<bool>("cracked", dict, patch.size());
    }

    if (dict.found("tractionN"))
    {
        tractionN_ = scalarField("tractionN", dict, patch.size());
    }

    if (dict.found("tractionS"))
    {
        tractionS_ = scalarField("tractionS", dict, patch.size());
    }

    if (dict.found("oldTractionN"))
    {
        oldTractionN_ = scalarField("oldTractionN", dict, patch.size());
    }

    if (dict.found("oldTractionS"))
    {
        oldTractionS_ = scalarField("oldTractionS", dict, patch.size());
    }

    if (dict.found("deltaN"))
    {
        deltaN_ = scalarField("deltaN", dict, patch.size());
    }

    if (dict.found("deltaS"))
    {
        deltaS_ = scalarField("deltaS", dict, patch.size());
    }

    if (dict.found("oldDeltaN"))
    {
        oldDeltaN_ = scalarField("oldDeltaN", dict, patch.size());
    }

    if (dict.found("oldDeltaS"))
    {
        oldDeltaS_ = scalarField("oldDeltaS", dict, patch.size());
    }

    if (dict.found("deltaEff"))
    {
        deltaEff_ = scalarField("deltaEff", dict, patch.size());
    }

    if (dict.found("unloadingDeltaEff"))
    {
        unloadingDeltaEff_ =
            scalarField("unloadingDeltaEff", dict, patch.size());
    }

    if (dict.found("GI"))
    {
        GI_ = scalarField("GI", dict, patch.size());
    }

    if (dict.found("GII"))
    {
        GII_ = scalarField("GII", dict, patch.size());
    }

    if (dict.found("oldGI"))
    {
        oldGI_ = scalarField("oldGI", dict, patch.size());
    }

    if (dict.found("oldGII"))
    {
        oldGII_ = scalarField("oldGII", dict, patch.size());
    }
}


// Construct as a copy
Foam::variableMixedModeCohesiveZoneModel::variableMixedModeCohesiveZoneModel
(
    const variableMixedModeCohesiveZoneModel& czm
)
:
    cohesiveZoneModel(czm),
    sigmaMax_(czm.sigmaMax_),
    tauMax_(czm.tauMax_),
    GIc_(czm.GIc_),
    GIIc_(czm.GIIc_),
    cracked_(czm.cracked_),
    tractionN_(czm.tractionN_),
    oldTractionN_(czm.oldTractionN_),
    tractionS_(czm.tractionS_),
    oldTractionS_(czm.oldTractionS_),
    deltaN_(czm.deltaN_),
    oldDeltaN_(czm.oldDeltaN_),
    deltaS_(czm.deltaS_),
    oldDeltaS_(czm.oldDeltaS_),
    deltaEff_(czm.deltaEff_),
    unloadingDeltaEff_(czm.unloadingDeltaEff_),
    GI_(czm.GI_),
    oldGI_(czm.oldGI_),
    GII_(czm.GII_),
    oldGII_(czm.oldGII_),
    penaltyScale_(czm.penaltyScale_),
    penaltyFactor_(czm.penaltyFactor_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::variableMixedModeCohesiveZoneModel::~variableMixedModeCohesiveZoneModel()
{}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::variableMixedModeCohesiveZoneModel::crackingAndDamage() const
{
    tmp<scalarField> tcrackingAndDamage
    (
     new scalarField(patch().size(), 1.0)
    );
    scalarField& cad = tcrackingAndDamage.ref();

    forAll(cad, faceI)
    {
        if (cracked_[faceI])
        {
            cad[faceI] = 2.0;
        }
    }

    return tcrackingAndDamage;
}


void Foam::variableMixedModeCohesiveZoneModel::autoMap
(
    const fvPatchFieldMapper& m
)
{
    {
        scalarField cracked(cracked_.size(), 1.0);
        forAll(cracked, i)
        {
            cracked[i] = cracked_[i] ? 1.0 : 0.0;
        }
        m(cracked, cracked);
        cracked_.resize(cracked.size());
        forAll(cracked, i)
        {
            cracked_[i] = cracked_[i] > 0.5;
        }
    }
    const label nNewFaces = cracked_.size() - tractionN_.size();

    m(tractionN_, tractionN_);
    m(oldTractionN_, oldTractionN_);

    m(tractionS_, tractionS_);
    m(oldTractionS_, oldTractionS_);

    m(deltaN_, deltaN_);
    m(oldDeltaN_, oldDeltaN_);

    m(deltaS_, deltaS_);
    m(oldDeltaS_, oldDeltaS_);

    m(deltaEff_, deltaEff_);
    m(unloadingDeltaEff_, unloadingDeltaEff_);

    m(GI_, GI_);
    m(oldGI_, oldGI_);

    m(GII_, GII_);
    m(oldGII_, oldGII_);


    // Only perform mapping if the number of faces on the patch has changed

    if (nNewFaces > 0 && isA<directFvPatchFieldMapper>(m))
    {
        // Reset values on new faces to zero
        // Note: the method below is used to find which faces are new on the
        // patch

        const directFvPatchFieldMapper& dm =
            dynamicCast<const directFvPatchFieldMapper&>(m);

        const labelList& addressing = dm.addressing();
        const label patchSize = patch().size();

        if (patchSize == 1 && nNewFaces == 1)
        {
            label i = 0;

            cracked_[i] = false;
            tractionN_[i] = 0.0;
            oldTractionN_[i] = 0.0;
            tractionS_[i] = 0.0;
            oldTractionS_[i] = 0.0;
            deltaN_[i] = 0.0;
            oldDeltaN_[i] = 0.0;
            deltaS_[i] = 0.0;
            oldDeltaS_[i] = 0.0;
            deltaEff_[i] = 0.0;
            unloadingDeltaEff_[i] = 0.0;
            GI_[i] = 0.0;
            oldGI_[i] = 0.0;
            GII_[i] = 0.0;
            oldGII_[i] = 0.0;
        }
        else if (patchSize == 2 && nNewFaces == 1)
        {
            label i = 1;

            cracked_[i] = false;
            tractionN_[i] = 0.0;
            oldTractionN_[i] = 0.0;
            tractionS_[i] = 0.0;
            oldTractionS_[i] = 0.0;
            deltaN_[i] = 0.0;
            oldDeltaN_[i] = 0.0;
            deltaS_[i] = 0.0;
            oldDeltaS_[i] = 0.0;
            deltaEff_[i] = 0.0;
            unloadingDeltaEff_[i] = 0.0;
            GI_[i] = 0.0;
            oldGI_[i] = 0.0;
            GII_[i] = 0.0;
            oldGII_[i] = 0.0;
        }
        else if (patchSize == 2 && nNewFaces == 2)
        {
            label i = 0;

            cracked_[i] = false;
            tractionN_[i] = 0.0;
            oldTractionN_[i] = 0.0;
            tractionS_[i] = 0.0;
            oldTractionS_[i] = 0.0;
            deltaN_[i] = 0.0;
            oldDeltaN_[i] = 0.0;
            deltaS_[i] = 0.0;
            oldDeltaS_[i] = 0.0;
            deltaEff_[i] = 0.0;
            unloadingDeltaEff_[i] = 0.0;
            GI_[i] = 0.0;
            oldGI_[i] = 0.0;
            GII_[i] = 0.0;
            oldGII_[i] = 0.0;

            i = 1;

            cracked_[i] = false;
            tractionN_[i] = 0.0;
            oldTractionN_[i] = 0.0;
            tractionS_[i] = 0.0;
            oldTractionS_[i] = 0.0;
            deltaN_[i] = 0.0;
            oldDeltaN_[i] = 0.0;
            deltaS_[i] = 0.0;
            oldDeltaS_[i] = 0.0;
            deltaEff_[i] = 0.0;
            unloadingDeltaEff_[i] = 0.0;
            GI_[i] = 0.0;
            oldGI_[i] = 0.0;
            GII_[i] = 0.0;
            oldGII_[i] = 0.0;
        }
        else
        {
            for (label i = 1; i < patchSize; i++)
            {
                if (addressing[i] == 0)
                {
                    cracked_[i] = false;
                    tractionN_[i] = 0.0;
                    oldTractionN_[i] = 0.0;
                    tractionS_[i] = 0.0;
                    oldTractionS_[i] = 0.0;
                    deltaN_[i] = 0.0;
                    oldDeltaN_[i] = 0.0;
                    deltaS_[i] = 0.0;
                    oldDeltaS_[i] = 0.0;
                    deltaEff_[i] = 0.0;
                    unloadingDeltaEff_[i] = 0.0;
                    GI_[i] = 0.0;
                    oldGI_[i] = 0.0;
                    GII_[i] = 0.0;
                    oldGII_[i] = 0.0;
                }
            }
        }
    }
}


void Foam::variableMixedModeCohesiveZoneModel::rmap
(
    const solidCohesiveFvPatchVectorField& sc,
    const labelList& addr
)
{
    const variableMixedModeCohesiveZoneModel& czm =
        refCast<const variableMixedModeCohesiveZoneModel>(sc.cohesiveZone());

    cracked_.rmap(czm.cracked_, addr);
    tractionN_.rmap(czm.tractionN_, addr);
    oldTractionN_.rmap(czm.oldTractionN_, addr);
    tractionS_.rmap(czm.tractionS_, addr);
    oldTractionS_.rmap(czm.oldTractionS_, addr);
    deltaN_.rmap(czm.deltaN_, addr);
    oldDeltaN_.rmap(czm.oldDeltaN_, addr);
    deltaS_.rmap(czm.deltaS_, addr);
    oldDeltaS_.rmap(czm.oldDeltaS_, addr);
    unloadingDeltaEff_.rmap(czm.unloadingDeltaEff_, addr);
    deltaEff_.rmap(czm.deltaEff_, addr);
    GI_.rmap(czm.GI_, addr);
    oldGI_.rmap(czm.oldGI_, addr);
    GII_.rmap(czm.GII_, addr);
    oldGII_.rmap(czm.oldGII_, addr);
}


void Foam::variableMixedModeCohesiveZoneModel::updateOldFields()
{
    // Update old tractions
    oldTractionN_ = tractionN_;
    oldTractionS_ = tractionS_;

    // Update energies
    oldGI_ = GI_;
    oldGII_ = GII_;

    // Update old deltas
    oldDeltaN_ = deltaN_;
    oldDeltaS_ = deltaS_;

    // Update unloading delta effective
    // It is the maximum of the  delta effective and the previous
    // unloading delta effective
    unloadingDeltaEff_ = max(unloadingDeltaEff_, deltaEff_);
}


void Foam::variableMixedModeCohesiveZoneModel::updateTraction
(
    vectorField& traction,
    const vectorField& delta
)
{
    // Unit normal vectors
    const vectorField& n = patch().patch().faceNormals();

    // Check if the penaltyFactor needs to be calculated
    if (penaltyFactor_ < SMALL)
    {
        calcPenaltyFactor();
    }

    // Note: deltaN, deltaS, tractionN, and tractionS have been updated by
    // updateEnergy, which is called before updateTraction

    const scalar faceGIc = GIc_.value();
    const scalar faceGIIc = GIIc_.value();

    const scalar faceSigmaMax = sigmaMax_.value();
    const scalar faceTauMax = tauMax_.value();

    // Set tractions for each face

    forAll(traction, faceI)
    {
        vector& faceTrac = traction[faceI];
        bool& faceCracked = cracked_[faceI];

        const vector& faceDelta = delta[faceI];
        const scalar faceDeltaN = deltaN_[faceI];
        const scalar faceDeltaS = deltaS_[faceI];
        const scalar faceDeltaEff = deltaEff_[faceI];
        const scalar faceUnloadingDeltaEff = unloadingDeltaEff_[faceI];

        const vector& faceN = n[faceI];

        const scalar faceGI = GI_[faceI];
        const scalar faceGII = GII_[faceI];

        // Check propagation criterion for new cracked faces
        if (!faceCracked && ((faceGI/faceGIc) + (faceGII/faceGIIc)) > 1.0)
        {
            Pout<< "Face " << faceI << " is fully cracked" << endl;
            faceCracked = true;
        }

        // The phase can be in one of three phases:
        //    1) cracked
        //    2) loading
        //    3) unloading

        // Calculate tractions
        if (faceCracked)
        {
            // Cracked faces may have a contact traction
            if (faceDeltaN < 0.0)
            {
                faceTrac = faceDeltaN*penaltyFactor_*faceN;
            }
            else
            {
                faceTrac = vector::zero;
            }
        }
        else if (faceDeltaEff >= faceUnloadingDeltaEff)
        {
            // Loading phase: i.e. further fracture energy is being dissipated.
            // If the effective delta is greater than unloadingDeltaEff then
            // there is loading

            if (faceDeltaN > 0.0)
            {
                // Face in loading damage phase which is not in contact

                // Set traction in a fixed point iteration manner to force
                // (tN/sigmaMax)^2 + (tS/tauMax)^2 = 1
                // while also allowing varying mode mix by assuming colinear
                // traction
                // and delta

                // Dugdale style CZM assumed

                // Note: previous method called
                // mechanical.cohLaw().damageTractions

                // Set traction in damage zone
                faceTrac =
                    damageTractionN
                    (
                        faceDeltaN, faceDeltaS, faceSigmaMax, faceTauMax, faceN
                    )
                  + damageTractionS
                    (
                        faceDeltaN, faceDeltaS, faceSigmaMax, faceTauMax, faceN,
                        faceDelta
                    );
            }
            else
            {
                // Face in loading damage phase which is in contact i.e. the
                // face may have just initiated and the contactOffset may be
                // greater than zero i.e. 0 < faceDeltaN < contactOffset

                // Set traction as the sum of normal contact and shear damage
                // Note: we set fricCoeff as 0.0 because the shear is set from
                // the damage shear traction
                faceTrac =
                    faceDeltaN*penaltyFactor_*faceN
                  + damageTractionS
                    (
                        faceDeltaN, faceDeltaS, faceSigmaMax, faceTauMax, faceN,
                        faceDelta
                    );
            }
        }
        else
        {
            // Unloading

            if (faceDeltaN < 0.0)
            {
                // Face is damaged and in contact and has yet to dissipate any
                // fracture energy i.e. the unloadingFaceDEff is zero
                faceTrac = faceDeltaN*penaltyFactor_*faceN;
            }
            else
            {
                // The  effective faceDelta is less than the unloading effective
                // faceDelta, and greater than the contact offset

                // We have two choice for unloading:
                // (a) ductile approach
                //         unload with initial stiffness (which is infinite)
                //         similar to ductilve metal stress-strain curve
                //         as we use infinite intial stiffness, this means to
                //         elastic energy is recovered and we immediately start
                //         loading in opposite direction
                // (b) brittle approach
                //         unload straight back to the origin
                //         this means we recover elastic energy
                //
                // Approach (b) is numerically "nicer"; however, for Dugdale
                // cohesive zone this implys that only half the energy is
                // dissipated just before failure and then the other half is
                // dissipated at failure. This may not be an issue in practice
                // but it does not really make sense.
                // Obviously it is fine for a linear cohesive zone as the energy
                // is smoothly dissipated up to failure. For now, we will
                // implement approach (b), but this requires more thinking...

                // Reduce traction linearly with the reduction in faceDelta
                if (faceUnloadingDeltaEff > SMALL)
                {
                    // Set tractions in damage zone and scale them due to
                    // unloading
                    faceTrac =
                        (faceDeltaEff/faceUnloadingDeltaEff)
                        *(
                            damageTractionN
                            (
                                faceDeltaN, faceDeltaS,
                                faceSigmaMax, faceTauMax,
                                faceN
                            )
                            + damageTractionS
                            (
                                faceDeltaN, faceDeltaS,
                                faceSigmaMax, faceTauMax,
                                faceN,
                                faceDelta
                            )
                        );
                }
                else
                {
                    faceTrac = vector::zero;
                }
            }
        }
    }

    // Note: we do not under-relax the traction field because instead the
    // faceDeltas are under-relaxed in solidCohesive
}

void Foam::variableMixedModeCohesiveZoneModel::updateEnergy
(
    const vectorField& traction,
    const vectorField& delta
)
{
    // Unit normal vectors
    const vectorField& n = patch().patch().faceNormals();

    // Update normal delta
    deltaN_ = n & delta;

    // Update shear delta
    deltaS_ = mag((I - sqr(n)) & delta);

    // Update delta effective, where only positive deltaN is considered
    deltaEff_ = Foam::sqrt(pow(max(deltaN_, 0.0), 2) + pow(deltaS_, 2));

    // Update normal traction
    tractionN_ = n & traction;

    // Update shear traction
    tractionS_ = mag((I - sqr(n)) & traction);


    // The dissipated fracture energy is calculated by integrating the product
    // the tractions and the deltas using the trapezoidal rule

    forAll(delta, faceI)
    {
        const scalar faceDeltaN = deltaN_[faceI];
        const scalar faceDeltaS = deltaS_[faceI];
        const scalar faceOldDeltaN = oldDeltaN_[faceI];
        const scalar faceOldDeltaS = oldDeltaS_[faceI];
        const scalar faceDeltaEff = deltaEff_[faceI];
        const scalar faceUnloadingDeltaEff = unloadingDeltaEff_[faceI];

        const scalar faceTracN = tractionN_[faceI];
        const scalar faceTracS = tractionS_[faceI];
        const scalar faceOldTracN = oldTractionN_[faceI];
        const scalar faceOldTracS = oldTractionS_[faceI];

        const bool faceCracked = cracked_[faceI];

        scalar& faceGI = GI_[faceI];
        scalar& faceGII = GII_[faceI];

        const scalar faceOldGI = oldGI_[faceI];
        const scalar faceOldGII = oldGII_[faceI];

        if (!faceCracked && (faceDeltaEff > faceUnloadingDeltaEff))
        {
            // If the average normal stress is tensile
            if ((faceTracN + faceOldTracN) > 0.0)
            {
                // Integrate using trapezoidal rule
                faceGI =
                    faceOldGI
                  + (
                        0.5*(faceTracN + faceOldTracN)
                       *(faceDeltaN - faceOldDeltaN)
                    );
            }
            else
            {
                // No modeI energy dissipated if the face is in compression
                faceGI = faceOldGI;
            }

            // Mode II - trapezoidal rule
            faceGII =
                faceOldGII
              + (
                    0.5*(faceTracS + faceOldTracS)*(faceDeltaS - faceOldDeltaS)
                );
        }
    }
}


Foam::tmp<Foam::surfaceScalarField>
Foam::variableMixedModeCohesiveZoneModel::initiationTractionFraction() const
{
    // Reference to the mesh
    const fvMesh& mesh = this->mesh();

    // Face unit normals
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // Lookup the traction field from the solver: this should be up-to-date
    const surfaceVectorField& traction = meshTraction();

    // Calculate normal traction
    const surfaceScalarField normalTrac
    (
        max(dimensionedScalar("zero", dimPressure, 0.0), n & traction)
    );

    // Calculate shear traction
    const surfaceScalarField shearTrac(mag((I - sqr(n)) & traction));

    // The tractionFraction >= 1.0 for faces to be broken and added to the
    // cohesive zone
    return
        tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "tractionFraction",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                sqrt(sqr(normalTrac/sigmaMax_) + sqr(shearTrac/tauMax_))
            )
        );
}


void Foam::variableMixedModeCohesiveZoneModel::write(Ostream& os) const
{
    writeEntry(os, "type", type());

    writeEntry(os, "sigmaMax", sigmaMax_.value());
    writeEntry(os, "tauMax", tauMax_.value());
    writeEntry(os, "GIc", GIc_.value());
    writeEntry(os, "GIIc", GIIc_.value());

    writeEntry(os, "cracked", cracked_);

    writeEntry(os, "tractionN", tractionN_);
    writeEntry(os, "tractionS", tractionS_);
    writeEntry(os, "oldTractionN", oldTractionN_);
    writeEntry(os, "oldTractionS", oldTractionS_);

    writeEntry(os, "deltaN", deltaN_);
    writeEntry(os, "deltaS", deltaS_);
    writeEntry(os, "oldDeltaN", oldDeltaN_);
    writeEntry(os, "oldDeltaS", oldDeltaS_);

    writeEntry(os, "deltaEff", deltaEff_);
    writeEntry(os, "unloadingDeltaEff", unloadingDeltaEff_);

    writeEntry(os, "GI", GI_);
    writeEntry(os, "GII", GII_);
    writeEntry(os, "oldGI", oldGI_);
    writeEntry(os, "oldGII", oldGII_);

    writeEntry(os, "penaltyScale", penaltyScale_);
    writeEntry(os, "penaltyFactor", penaltyFactor_);
}


// ************************************************************************* //
