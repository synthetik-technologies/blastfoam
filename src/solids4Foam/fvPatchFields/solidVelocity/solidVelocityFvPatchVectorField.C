/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "solidVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

void solidVelocityFvPatchVectorField::makeInterp() const
{
    if (interpPtr_.valid())
    {
        FatalErrorIn
        (
            "void solidVelocityFvPatchVectorField::makeInterp() const"
        ) << "pointer already set" << abort(FatalError);
    }

    interpPtr_.set(new primitivePatchInterpolation(patch().patch()));
}


primitivePatchInterpolation& solidVelocityFvPatchVectorField::interp()
{
    if (interpPtr_.empty())
    {
        makeInterp();
    }

    return interpPtr_();
}


void solidVelocityFvPatchVectorField::setPointDisplacement
(
    const vectorField& faceDisp
)
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if
    (
        mesh.foundObject<pointVectorField>
        (
            "point" + internalField().name()
        )
    )
    {
        const pointVectorField& pointD =
            mesh.lookupObject<pointVectorField>
            (
                "point" + internalField().name()
            );

        // Check if the boundary is fixedValue
        if
        (
            pointD.boundaryField()[patch().index()].type()
         == fixedValuePointPatchVectorField::typeName
        )
        {
            // Use const_cast to set boundary condition
            fixedValuePointPatchVectorField& patchPointD =
                refCast<fixedValuePointPatchVectorField>
                (
                    const_cast<pointVectorField&>
                    (
                        pointD
                    ).boundaryFieldRef()[patch().index()]
                );

            // Interpolate face values to the points
            patchPointD == interp().faceToPointInterpolate(faceDisp);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidVelocityFvPatchVectorField::solidVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    velocity_(p.size(), vector::zero),
    velocitySeries_(),
    interpPtr_(NULL)
{}


solidVelocityFvPatchVectorField::solidVelocityFvPatchVectorField
(
    const solidVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    velocity_(mapper(ptf.velocity_)),
    velocitySeries_(ptf.velocitySeries_, false),
    interpPtr_(NULL)
{}


solidVelocityFvPatchVectorField::solidVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    velocity_(p.size(), vector::zero),
    velocitySeries_(),
    interpPtr_(NULL)
{
    Info<< "Creating " << type() << " boundary condition" << endl;

    // Read velocity
    if (dict.found("velocity"))
    {
        velocity_ = vectorField("velocity", dict, p.size());
    }
    else if (dict.found("velocitySeries"))
    {
        Info<< "    velocity is time-varying" << endl;
        velocitySeries_ = Function1<vector>::New("velocitySeries", dict);

        fvPatchField<vector>::operator==
        (
            velocitySeries_->value(this->db().time().timeOutputValue())
        );
    }
    else
    {
        FatalErrorIn(type() + "::solidVelocityFvPatchVectorField(...)")
            << "Either 'velocity' or 'velocitySeries' should be specified!"
            << abort(FatalError);
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }
}


solidVelocityFvPatchVectorField::solidVelocityFvPatchVectorField
(
    const solidVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    velocity_(pivpvf.velocity_),
    velocitySeries_(pivpvf.velocitySeries_, false),
    interpPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void solidVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);

    m(velocity_, velocity_);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidVelocityFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const solidVelocityFvPatchVectorField& dmptf =
       refCast<const solidVelocityFvPatchVectorField>(ptf);

    velocity_.rmap(dmptf.velocity_, addr);
}


void solidVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Check if the velocity is time-varying
    if (velocitySeries_.valid())
    {
        velocity_ = velocitySeries_->value(this->db().time().timeOutputValue());
    }

    vectorField disp = vectorField(patch().size(), vector::zero);

    if (internalField().name() == "DD")
    {
        // Incremental approach, so we wil set the increment of displacement for
        // this time-step
        disp = velocity_*db().time().deltaTValue();
    }
    else
    {
        // Lookup the old time total displacement
        const volVectorField& Dold =
            db().lookupObject<volVectorField>("D").oldTime();

        // The new total displacement is equal to Dold plus the increment of
        // displacement based on the current velocity and time-step
        disp =
            Dold.boundaryField()[patch().index()]
          + velocity_*db().time().deltaTValue();
    }

    // Set the displacement (or displacement increment) on the patch
    fvPatchField<vector>::operator==(disp);
    fixedValueFvPatchVectorField::updateCoeffs();

    // If the corresponding point displacement field has a fixedValue type
    // boundary condition, then we wil update it
    setPointDisplacement(disp);
}


Foam::tmp<Foam::Field<vector> > solidVelocityFvPatchVectorField::snGrad() const
{
    // fixedValue snGrad with no correction
    // return (*this - patchInternalField())*this->patch().deltaCoeffs();

    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + internalField().name() + ")"
        );

    // Unit normals
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Correction vectors
    const vectorField k(delta - n*(n&delta));

    return
    (
        *this - (patchInternalField() + (k & gradField.patchInternalField()))
    )*patch().deltaCoeffs();
}

tmp<Field<vector> >
solidVelocityFvPatchVectorField::gradientBoundaryCoeffs() const
{
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + internalField().name() + ")"
        );

    vectorField n(this->patch().nf());
    vectorField delta(this->patch().delta());

    //- correction vector
    vectorField k(delta - n*(n&delta));

    return
    (
        this->patch().deltaCoeffs()
       *(*this - (k & gradField.patchInternalField()))
    );
}

void solidVelocityFvPatchVectorField::write(Ostream& os) const
{
    if (velocitySeries_.valid())
    {
        writeEntry(os, "velocitySeries", velocitySeries_());
    }
    else
    {
        writeEntry(os, "velocity", velocity_);
    }

    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    solidVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
