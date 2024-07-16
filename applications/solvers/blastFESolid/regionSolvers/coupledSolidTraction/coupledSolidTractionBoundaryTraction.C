
#include "coupledSolidTractionBoundaryTraction.H"
#include "pointFields.H"
#include "volFields.H"
#include "compressibleMomentumTransportModel.H"
#include "incompressibleMomentumTransportModel.H"
#include "globalPolyBoundaryMesh.H"
#include "coupledGlobalPolyPatch.H"
#include "vtkWritePolyData.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace boundaryTractions
{
    defineTypeNameAndDebug(coupledSolidTraction, 0);
    addToRunTimeSelectionTable
    (
        boundaryTraction,
        coupledSolidTraction,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::symmTensorField>
Foam::boundaryTractions::coupledSolidTraction::viscousStress
(
    const polyMesh& mesh,
    const polyPatch& patch
) const
{
    typedef compressibleMomentumTransportModel cmpTurbModel;
    typedef incompressibleMomentumTransportModel icoTurbModel;

    if (mesh.foundObject<volSymmTensorField>("devTau"))
    {
        return
            mesh.lookupObject<volSymmTensorField>
            (
                "devTau"
            ).boundaryField()[patch.index()];
    }
    else if (mesh.foundObject<cmpTurbModel>(cmpTurbModel::typeName))
    {
        const cmpTurbModel& turb
        (
            mesh.lookupObject<cmpTurbModel>(cmpTurbModel::typeName)
        );

        return turb.devTau()().boundaryField()[patch.index()];
    }
    else if (mesh.foundObject<volSymmTensorField>("devSigma"))
    {

        return
            mesh.lookupObject<volSymmTensorField>
            (
                "devSigma"
            ).boundaryField()[patch.index()]
           *rho(mesh, patch);

    }
    else if (mesh.foundObject<icoTurbModel>(icoTurbModel::typeName))
    {
        const icoTurbModel& turb
        (
            mesh.lookupObject<icoTurbModel>(icoTurbModel::typeName)
        );

        return
            turb.devSigma()().boundaryField()[patch.index()]
           *rho(mesh, patch);
    }
    else
    {
        // NotImplemented;
        // For laminar flows get the velocity
        // const fvPatchVectorField& Up
        // (
        //     patch.lookupPatchField<volVectorField, vector>("U")
        // );

        // return mu(mesh, patch)*Up.snGrad();
        return tmp<symmTensorField>
        (
            new symmTensorField(patch.size(), Zero)
        );
    }
}


Foam::tmp<Foam::scalarField>
Foam::boundaryTractions::coupledSolidTraction::rho
(
    const polyMesh& mesh,
    const polyPatch& patch
) const
{
    if (mesh.foundObject<volScalarField>("rho"))
    {
        return
            mesh.lookupObject<volScalarField>
            (
                "rho"
            ).boundaryField()[patch.index()];
    }
    if (mesh.foundObject<volScalarField>("thermo:rho"))
    {
        return
            mesh.lookupObject<volScalarField>
            (
                "thermo::rho"
            ).boundaryField()[patch.index()];
    }
    else if (mesh.foundObject<IOdictionary>("transportProperties"))
    {
        const IOdictionary& transportProperties =
            mesh.lookupObject<IOdictionary>("transportProperties");

        return tmp<scalarField>
        (
            new scalarField
            (
                patch.size(),
                transportProperties.lookup<scalar>("rho")
            )
        );
    }
    else
    {
        return tmp<scalarField>(new scalarField(patch.size(), 1.0));
    }
}


Foam::tmp<Foam::scalarField>
Foam::boundaryTractions::coupledSolidTraction::mu
(
    const polyMesh& mesh,
    const polyPatch& patch
) const
{
    if (mesh.foundObject<volScalarField>("thermo:mu"))
    {
        return
            mesh.lookupObject<volScalarField>
            (
                "thermo::mu"
            ).boundaryField()[patch.index()];
    }
    else if (mesh.foundObject<volScalarField>("mu"))
    {
        return
            mesh.lookupObject<volScalarField>
            (
                "mu"
            ).boundaryField()[patch.index()];
    }
    else if (mesh.foundObject<IOdictionary>("transportProperties"))
    {
        const IOdictionary& transportProperties =
            mesh.lookupObject<IOdictionary>("transportProperties");

        if (transportProperties.found("nu"))
        {
             return
                rho(mesh, patch)
               *transportProperties.lookup<scalar>("nu");
        }

        return tmp<scalarField>
        (
            new scalarField
            (
                patch.size(),
                transportProperties.lookup<scalar>("mu")
            )
        );
    }
    else
    {
        return tmp<scalarField>(new scalarField(patch.size(), 0.0));
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryTractions::coupledSolidTraction::coupledSolidTraction
(
    const word& name,
    const dictionary& dict,
    const feMesh1& femesh,
    const pointVectorField* DPtr
)
:
    boundaryTraction(name, dict, femesh, DPtr),
    pName_(dict.lookupOrDefault("pName", word("p"))),
    pRef_(dict.lookup<scalar>("pRef")),
    cutOffPressure_(dict.lookupOrDefault("cutOffPressure", 0.0))
{}


Foam::boundaryTractions::coupledSolidTraction::coupledSolidTraction
(
    const coupledSolidTraction& cstbt
)
:
    boundaryTraction(cstbt),
    pName_(cstbt.pName_),
    pRef_(cstbt.pRef_),
    cutOffPressure_(cstbt.cutOffPressure_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boundaryTractions::coupledSolidTraction::~coupledSolidTraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boundaryTractions::coupledSolidTraction::update()
{
    // Get the coupling information from the mappedPatchBase
    const coupledGlobalPolyPatch& cgpp =
        globalPolyBoundaryMesh::New(mesh_)
        (
            mesh_.boundaryMesh()[patchID_]
        );
    const polyMesh& nbrMesh = cgpp.sampleMesh();
    const coupledGlobalPolyPatch& samplePatch = cgpp.samplePatch();
    const coupledGlobalPolyPatch& nbrCgpp =
        globalPolyBoundaryMesh::New(nbrMesh)(cgpp.samplePatch().patch());
    const polyPatch& samplePolyPatch = nbrCgpp.patch();
    const label samplePatchi = samplePolyPatch.index();

    //- Lookup pressure fields
    const volScalarField& pNbr =
        nbrMesh.lookupObject<volScalarField>(pName_);
    scalarField ppNbr(pNbr.boundaryField()[samplePatchi]);

    if (pNbr.dimensions() != dimPressure)
    {
        ppNbr *= rho(nbrMesh, samplePolyPatch);
    }

    if (cutOffPressure_ > 0)
    {
        forAll(ppNbr, i)
        {
            if (mag(ppNbr[i]) < cutOffPressure_)
            {
                ppNbr[i] = 0.0;
            }
        }
    }

//     if (debug > 1 || (debug && mesh_.time().outputTime()))
//     {
//         Field<scalar> pfGlobal(samplePatch.patchFaceToGlobal(ppNbr));
//         Field<scalar> pfInterp
//         (
//             cgpp.patchToPatchInterpolator().transferFaces
//             (
//                 samplePatch.globalPatch(),
//                 pfGlobal
//             )
//         );
//
//         if (Pstream::master())
//         {
//             fileName path
//             (
//                 mesh_.time().globalPath()
//                /"VTK"
//                /mesh_.time().timeName()
//             );
//             mkDir(path);
//             vtkWritePolyData::write
//             (
//                 path/"p_interpolated.vtk",
//                 "p",
//                 true,
//                 cgpp.globalPatch().points(),
//                 labelList(),
//                 edgeList(),
//                 cgpp.globalPatch(),
//                 "p",
//                 false,
//                 pfInterp
//
//             );
//             vtkWritePolyData::write
//             (
//                 path/"p_actual.vtk",
//                 "p",
//                 true,
//                 samplePatch.globalPatch().points(),
//                 labelList(),
//                 edgeList(),
//                 samplePatch.globalPatch(),
//                 "p",
//                 false,
//                 pfGlobal
//             );
//         }
//     }

    pressure_ = cgpp.faceToPoint(nbrCgpp.faceInterpolate(ppNbr));

    // Map viscous stress
    sigma_ = cgpp.faceToPoint
    (
        nbrCgpp.faceInterpolate(viscousStress(nbrMesh, samplePolyPatch))
    );
}


Foam::scalar Foam::boundaryTractions::coupledSolidTraction::pressure
(
    const labelList& nodeLabels,
    const scalarList& shape,
    const vector& n
) const
{
    scalar p = 0.0;
    forAll(shape, i)
    {
        p += shape[i]*pressure_[nodeLabels[i]];
    }
    return p - pRef_;
}


Foam::vector Foam::boundaryTractions::coupledSolidTraction::traction
(
    const labelList& nodeLabels,
    const scalarList& shape,
    const vector& n
) const
{
    symmTensor sigma(Zero);
    forAll(shape, i)
    {
        sigma += shape[i]*sigma_[nodeLabels[i]];
    }
    return sigma & n;
}


// ************************************************************************* //
