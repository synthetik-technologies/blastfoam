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

#include "modeICohesiveZoneModel.H"
#include "addToRunTimeSelectionTable.H"
#include "solidCohesiveFvPatchVectorField.H"
#include "directFvPatchFieldMapper.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(modeICohesiveZoneModel, 0)
    addToRunTimeSelectionTable
    (
        cohesiveZoneModel, modeICohesiveZoneModel, dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::vector Foam::modeICohesiveZoneModel::damageTractionN
(
    const scalar faceDeltaN,
    const scalar faceSigmaMax,
    const vector& faceN
) const
{
    return faceN*faceSigmaMax;
}


void Foam::modeICohesiveZoneModel::calcPenaltyFactor() const
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
Foam::modeICohesiveZoneModel::modeICohesiveZoneModel
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict
)
:
    cohesiveZoneModel(name, patch, dict),
    sigmaMax_("sigmaMax", dimPressure, dict),
    GIc_("GIc", dimensionSet(1, 0, -2, 0, 0, 0, 0), dict),
    cracked_(patch.size(), false),
    tractionN_(patch.size(), 0.0),
    oldTractionN_(patch.size(), 0.0),
    deltaN_(patch.size(), 0.0),
    oldDeltaN_(patch.size(), 0.0),
    deltaEff_(patch.size(), 0.0),
    unloadingDeltaEff_(patch.size(), 0.0),
    GI_(patch.size(), 0.0),
    oldGI_(patch.size(), 0.0),
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

    if (dict.found("oldTractionN"))
    {
        oldTractionN_ = scalarField("oldTractionN", dict, patch.size());
    }

    if (dict.found("deltaN"))
    {
        deltaN_ = scalarField("deltaN", dict, patch.size());
    }

    if (dict.found("oldDeltaN"))
    {
        oldDeltaN_ = scalarField("oldDeltaN", dict, patch.size());
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

    if (dict.found("oldGI"))
    {
        oldGI_ = scalarField("oldGI", dict, patch.size());
    }
}


// Construct as a copy
Foam::modeICohesiveZoneModel::modeICohesiveZoneModel
(
    const modeICohesiveZoneModel& czm
)
:
    cohesiveZoneModel(czm),
    sigmaMax_(czm.sigmaMax_),
    GIc_(czm.GIc_),
    cracked_(czm.cracked_),
    tractionN_(czm.tractionN_),
    oldTractionN_(czm.oldTractionN_),
    deltaN_(czm.deltaN_),
    oldDeltaN_(czm.oldDeltaN_),
    deltaEff_(czm.deltaEff_),
    unloadingDeltaEff_(czm.unloadingDeltaEff_),
    GI_(czm.GI_),
    oldGI_(czm.oldGI_),
    penaltyScale_(czm.penaltyScale_),
    penaltyFactor_(czm.penaltyFactor_)
{}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::modeICohesiveZoneModel::~modeICohesiveZoneModel()
{}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::modeICohesiveZoneModel::crackingAndDamage() const
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


void Foam::modeICohesiveZoneModel::autoMap(const fvPatchFieldMapper& m)
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
    m(deltaN_, deltaN_);
    m(oldDeltaN_, oldDeltaN_);
    m(deltaEff_, deltaEff_);
    m(unloadingDeltaEff_, unloadingDeltaEff_);
    m(GI_, GI_);
    m(oldGI_, oldGI_);

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
            deltaN_[i] = 0.0;
            oldDeltaN_[i] = 0.0;
            deltaEff_[i] = 0.0;
            unloadingDeltaEff_[i] = 0.0;
            GI_[i] = 0.0;
            oldGI_[i] = 0.0;
        }
        else if (patchSize == 2 && nNewFaces == 1)
        {
            label i = 1;

            cracked_[i] = false;
            tractionN_[i] = 0.0;
            oldTractionN_[i] = 0.0;
            deltaN_[i] = 0.0;
            oldDeltaN_[i] = 0.0;
            deltaEff_[i] = 0.0;
            unloadingDeltaEff_[i] = 0.0;
            GI_[i] = 0.0;
            oldGI_[i] = 0.0;
        }
        else if (patchSize == 2 && nNewFaces == 2)
        {
            label i = 0;

            cracked_[i] = false;
            tractionN_[i] = 0.0;
            oldTractionN_[i] = 0.0;
            deltaN_[i] = 0.0;
            oldDeltaN_[i] = 0.0;
            deltaEff_[i] = 0.0;
            unloadingDeltaEff_[i] = 0.0;
            GI_[i] = 0.0;
            oldGI_[i] = 0.0;

            i = 1;

            cracked_[i] = false;
            tractionN_[i] = 0.0;
            oldTractionN_[i] = 0.0;
            deltaN_[i] = 0.0;
            oldDeltaN_[i] = 0.0;
            deltaEff_[i] = 0.0;
            unloadingDeltaEff_[i] = 0.0;
            GI_[i] = 0.0;
            oldGI_[i] = 0.0;
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
                    deltaN_[i] = 0.0;
                    oldDeltaN_[i] = 0.0;
                    deltaEff_[i] = 0.0;
                    unloadingDeltaEff_[i] = 0.0;
                    GI_[i] = 0.0;
                    oldGI_[i] = 0.0;
                }
            }
        }
    }
}


void Foam::modeICohesiveZoneModel::rmap
(
    const solidCohesiveFvPatchVectorField& sc,
    const labelList& addr
)
{
    const modeICohesiveZoneModel& czm =
        refCast<const modeICohesiveZoneModel>(sc.cohesiveZone());

    cracked_.rmap(czm.cracked_, addr);
    tractionN_.rmap(czm.tractionN_, addr);
    oldTractionN_.rmap(czm.oldTractionN_, addr);
    deltaN_.rmap(czm.deltaN_, addr);
    oldDeltaN_.rmap(czm.oldDeltaN_, addr);
    unloadingDeltaEff_.rmap(czm.unloadingDeltaEff_, addr);
    deltaEff_.rmap(czm.deltaEff_, addr);
    GI_.rmap(czm.GI_, addr);
    oldGI_.rmap(czm.oldGI_, addr);
}


void Foam::modeICohesiveZoneModel::updateOldFields()
{
    // Update old tractions
    oldTractionN_ = tractionN_;

    // Update energies
    oldGI_ = GI_;

    // Update old deltas
    oldDeltaN_ = deltaN_;

    // Update unloading delta effective
    // It is the maximum of the  delta effective and the previous
    // unloading delta effective
    unloadingDeltaEff_ = max(unloadingDeltaEff_, deltaEff_);
}


void Foam::modeICohesiveZoneModel::updateTraction
(
    vectorField& traction,
    const vectorField& delta
)
{
    // Unit normal vectors
    const vectorField& n = patch().patch().faceNormals();

    // Note: deltaN, deltaS, tractionN, and tractionS have been updated by
    // updateEnergy, which is called before updateTraction

    const scalar faceGIc = GIc_.value();
    const scalar faceSigmaMax = sigmaMax_.value();

    // Check if the penaltyFactor needs to be calculated
    if (penaltyFactor_ < SMALL)
    {
        calcPenaltyFactor();
    }

    forAll(traction, faceI)
    {
        vector& faceTrac = traction[faceI];
        bool& faceCracked = cracked_[faceI];

        const scalar faceDeltaN = deltaN_[faceI];
        const scalar faceDeltaEff = deltaEff_[faceI];
        const scalar faceUnloadingDeltaEff = unloadingDeltaEff_[faceI];
        const vector& faceN = n[faceI];
        const scalar faceGI = GI_[faceI];

        // Check propagation criterion for new cracked faces
        if (!faceCracked && (faceGI/faceGIc) > 1.0)
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
            // Check for contact i.e. faces crossing over
            if (faceDeltaN < 0.0)
            {
                // Cracked faces may have a contact traction; we will use a
                // penalty approach and neglect friction
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
            // there is loading. Face in loading damage phase

            // Set traction in damage zone
            faceTrac = damageTractionN(faceDeltaN, faceSigmaMax, faceN);
        }
        else
        {
            // Unloading

            if (faceDeltaN < 0.0)
            {
                // Cracked faces may have a contact traction; we will use a
                // penalty approach and neglect friction
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
                        *damageTractionN(faceDeltaN, faceSigmaMax, faceN);
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


void Foam::modeICohesiveZoneModel::updateEnergy
(
    const vectorField& traction,
    const vectorField& delta
)
{
    // Unit normal vectors
    const vectorField& n = patch().patch().faceNormals();

    // Update normal delta
    deltaN_ = n & delta;

    // Update delta effective, where only positive deltaN is considered
    //deltaEff_ = Foam::sqrt(pow(max(deltaN_, 0.0), 2) + pow(deltaS_, 2));
    deltaEff_ = max(deltaN_, 0.0);

    // Update normal traction
    tractionN_ = n & traction;

    // The dissipated fracture energy is calculated by integrating the product
    // the tractions and the deltas using the trapezoidal rule
    forAll(delta, faceI)
    {
        const scalar faceDeltaN = deltaN_[faceI];
        const scalar faceOldDeltaN = oldDeltaN_[faceI];
        const scalar faceDeltaEff = deltaEff_[faceI];
        const scalar faceUnloadingDeltaEff = unloadingDeltaEff_[faceI];

        const scalar faceTracN = tractionN_[faceI];
        const scalar faceOldTracN = oldTractionN_[faceI];
        const bool faceCracked = cracked_[faceI];

        scalar& faceGI = GI_[faceI];
        const scalar faceOldGI = oldGI_[faceI];

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
                // No mode-I energy dissipated if the face is in compression
                faceGI = faceOldGI;
            }
        }
    }
}


Foam::tmp<Foam::surfaceScalarField>
Foam::modeICohesiveZoneModel::initiationTractionFraction() const
{
    // Reference to the mesh
    const fvMesh& mesh = this->mesh();

    // Face unit normals
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // Calculate the current mesh normal traction field
    const surfaceScalarField normalTraction(n & meshTraction());

    // Calculate the tractionFraction = traction/sigmaMax, where only tensile
    // normal tractions are considered.
    // The tractionFraction >= 1.0 for internal faces to be broken and added to
    // the cohesive zone
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
                 max
                (
                    dimensionedScalar("zero", dimPressure, 0.0),
                    normalTraction
                )/sigmaMax_
            )
        );
}


void Foam::modeICohesiveZoneModel::write(Ostream& os) const
{
    writeEntry(os, "type", type());

    writeEntry(os, "sigmaMax", sigmaMax_.value());
    writeEntry(os, "GIc", GIc_.value());

    writeEntry(os, "cracked", cracked_);
    writeEntry(os, "tractionN", tractionN_);
    writeEntry(os, "oldTractionN", oldTractionN_);
    writeEntry(os, "deltaN", deltaN_);
    writeEntry(os, "oldDeltaN", oldDeltaN_);
    writeEntry(os, "deltaEff", deltaEff_);
    writeEntry(os, "unloadingDeltaEff", unloadingDeltaEff_);
    writeEntry(os, "GI", GI_);
    writeEntry(os, "oldGI", oldGI_);

    writeEntry(os, "penaltyScale", penaltyScale_);
    writeEntry(os, "penaltyFactor", penaltyFactor_);
}


// ************************************************************************* //
