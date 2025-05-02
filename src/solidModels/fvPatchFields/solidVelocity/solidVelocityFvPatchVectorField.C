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
#include "valuePointPatchFields.H"


// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

void Foam::solidVelocityFvPatchVectorField::makeInterp() const
{
    if (interpPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    interpPtr_.set(new primitivePatchInterpolation(patch().patch()));
}


const Foam::solidVelocityFvPatchVectorField::primitivePatchInterpolation&
Foam::solidVelocityFvPatchVectorField::interp() const
{
    if (interpPtr_.empty())
    {
        makeInterp();
    }

    return interpPtr_();
}


void Foam::solidVelocityFvPatchVectorField::setPointDisplacement
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
            isA<valuePointPatchVectorField>
            (pointD.boundaryField()[patch().index()])
        )
        {
            // Use const_cast to set boundary condition
            valuePointPatchVectorField& patchPointD =
                refCast<valuePointPatchVectorField>
                (
                    const_cast<pointVectorField&>
                    (
                        pointD
                    ).boundaryFieldRef()[patch().index()]
                );
            Field<vector>& pD(patchPointD);

            // Interpolate face values to the points
            pD = interp().faceToPointInterpolate(faceDisp);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidVelocityFvPatchVectorField::solidVelocityFvPatchVectorField
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


Foam::solidVelocityFvPatchVectorField::solidVelocityFvPatchVectorField
(
    const solidVelocityFvPatchVectorField& svpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(svpvf, p, iF, mapper),
    velocity_(mapper(svpvf.velocity_)),
    velocitySeries_(svpvf.velocitySeries_, false),
    interpPtr_(NULL)
{}


Foam::solidVelocityFvPatchVectorField::solidVelocityFvPatchVectorField
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
    DebugInfo
        << "Creating " << type() << " boundary condition" << endl;

    // Read velocity
    if (dict.found("velocity"))
    {
        velocity_ = vectorField("velocity", dict, p.size());
    }
    else if (dict.found("velocitySeries"))
    {
        DebugInfo<< "    velocity is time-varying" << endl;
        velocitySeries_ = Function1<vector>::New
        (
            "velocitySeries",
            this->db().time().userUnits(),
            dimVelocity,
            dict
        );

        fvPatchField<vector>::operator==
        (
            velocitySeries_->value(this->db().time().value())
        );
    }
    else
    {
        FatalErrorInFunction
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


Foam::solidVelocityFvPatchVectorField::solidVelocityFvPatchVectorField
(
    const solidVelocityFvPatchVectorField& svpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(svpvf, iF),
    velocity_(svpvf.velocity_),
    velocitySeries_(svpvf.velocitySeries_, false),
    interpPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidVelocityFvPatchVectorField::map
(
    const fvPatchField<vector>& ptf,
    const fieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::map(ptf, mapper);

    const solidVelocityFvPatchVectorField& svpvf =
       refCast<const solidVelocityFvPatchVectorField>(ptf);

    mapper(velocity_, svpvf.velocity_);
}


void Foam::solidVelocityFvPatchVectorField::reset
(
    const fvPatchField<vector>& ptf
)
{
    fixedValueFvPatchVectorField::reset(ptf);

    const solidVelocityFvPatchVectorField& svpvf =
       refCast<const solidVelocityFvPatchVectorField>(ptf);

    velocity_.reset(svpvf.velocity_);
}


void Foam::solidVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Check if the velocity is time-varying
    if (velocitySeries_.valid())
    {
        velocity_ = velocitySeries_->value(this->db().time().value());
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


Foam::tmp<Foam::Field<Foam::vector>>
Foam::solidVelocityFvPatchVectorField::snGrad() const
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

Foam::tmp<Foam::Field<Foam::vector>>
Foam::solidVelocityFvPatchVectorField::gradientBoundaryCoeffs() const
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


void Foam::solidVelocityFvPatchVectorField::write(Ostream& os) const
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

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        solidVelocityFvPatchVectorField
    );
}


// ************************************************************************* //
