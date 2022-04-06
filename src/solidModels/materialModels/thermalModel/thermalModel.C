/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

Class
    thermalModel

Description
    Thermal  material properties for solids.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "thermalModel.H"
#include "fluidBlastThermo.H"
#include "solidBlastThermo.H"
#include "volFields.H"
#include "fvc.H"
#include "solidSubMeshes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermalModel, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalModel::thermalModel(const fvMesh& mesh, const bool isSolid)
:
    mesh_(mesh)
{
    IOdictionary thermophysicalProperties
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    if (mesh.foundObject<solidSubMeshes>(solidSubMeshes::typeName))
    {
        subsetMeshes_.set
        (
            &mesh.lookupObjectRef<solidSubMeshes>(solidSubMeshes::typeName)
        );
        const solidSubMeshes::meshSubsetList& subMeshes =
            subsetMeshes_->subMeshes();
        thermoModels_.setSize(subMeshes.size());

        if (isSolid)
        {
            thermoPtr_.set
            (
                solidBlastThermo::New
                (
                    mesh,
                    thermophysicalProperties.subDict
                    (
                        subMeshes[0].subMesh().name()
                    )
                ).ptr()
            );
        }
        else
        {
            thermoPtr_.set
            (
                fluidBlastThermo::New
                (
                    1,
                    mesh,
                    thermophysicalProperties.subDict
                    (
                        subMeshes[0].subMesh().name()
                    )
                ).ptr()
            );
        }
        forAll(subMeshes, i)
        {
            {
                volScalarField T
                (
                    IOobject
                    (
                        "T",
                        subMeshes[i].subMesh().time().timeName(),
                        subMeshes[i].subMesh()
                    ),
                    subMeshes[i].interpolate(thermoPtr_->T()),
                    calculatedFvPatchScalarField::typeName
                );
                T.write();
            }

            if (isSolid)
            {
                thermoModels_.set
                (
                    i,
                    solidBlastThermo::New
                    (
                        subMeshes[i].subMesh(),
                        thermophysicalProperties.subDict
                        (
                            subMeshes[i].subMesh().name()
                        )
                    ).ptr()
                );
            }
            else
            {
                thermoModels_.set
                (
                    i,
                    fluidBlastThermo::New
                    (
                        1,
                        subMeshes[i].subMesh(),
                        thermophysicalProperties.subDict
                        (
                            subMeshes[i].subMesh().name()
                        )
                    ).ptr()
                );
            }
        }
        correct();
    }
    else if (isSolid)
    {
        thermoPtr_.set
        (
            solidBlastThermo::New
            (
                mesh,
                thermophysicalProperties
            ).ptr()
        );
    }
    else
    {
        thermoPtr_.set
        (
            fluidBlastThermo::New
            (
                1,
                mesh,
                thermophysicalProperties
            ).ptr()
        );
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalModel::~thermalModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::thermalModel::rho() const
{
    return thermoPtr_->rho();
}


Foam::volScalarField& Foam::thermalModel::rho()
{
    return thermoPtr_->rho();
}


Foam::tmp<Foam::volScalarField> Foam::thermalModel::C() const
{
    return thermoPtr_->Cv();
}


Foam::tmp<Foam::volScalarField> Foam::thermalModel::k() const
{
    return thermoPtr_->kappa();
}


void Foam::thermalModel::correct()
{
    if (subsetMeshes_.valid())
    {
        // Accumulated subMesh fields and then map to the base mesh
        const solidSubMeshes::meshSubsetList& subMeshes =
            subsetMeshes_->subMeshes();
        PtrList<volScalarField> hes(subMeshes.size());
        PtrList<volScalarField> rhos(subMeshes.size());
        forAll(subMeshes, thermoi)
        {
            rhos.set
            (
                thermoi,
                volScalarField::New("rho", thermoModels_[thermoi].rho())
            );
            hes.set
            (
                thermoi,
                volScalarField::New("he", thermoModels_[thermoi].he())
            );
        }

        // Map subMesh fields to the base mesh
        subsetMeshes_->mapSubMeshVolFields<scalar>(rhos, thermoPtr_->rho());
        subsetMeshes_->mapSubMeshVolFields<scalar>(hes, thermoPtr_->he());

        // Clear subMesh fields
        rhos.clear();
        hes.clear();
    }

    thermoPtr_->correct();
    thermoPtr_->T().correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
