#include "fvCFD.H"
#include "pointFields.H"
#include "valuePointPatchFields.H"
#include "dynamicBlastFvMesh.H"
#include "extendedNLevelGlobalCellToCellStencils.H"
#include "Random.H"
#include "volPointInterpolation.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicBlastFvMesh.H"

    runTime.setDeltaT(1);
    runTime.setWriteInterval(1);

    volScalarField e
    (
        IOobject
        (
            "e",
            runTime.timeName(),
            mesh
        ),
        mesh,
        10.0
    );

    pointVectorField& pointDisplacement =
        mesh.lookupObjectRef<pointVectorField>("pointDisplacement");
    const volVectorField& cellDisplacement =
        mesh.lookupObject<volVectorField>("cellDisplacement");

    runTime++;


    forAll(pointDisplacement.boundaryField(), patchi)
    {
        pointPatchField<vector>& ppd =
            pointDisplacement.boundaryFieldRef()[patchi];
        if (isA<valuePointPatchField<vector>>(ppd))
        {
            valuePointPatchField<vector>& pvpd =
                dynamicCast<valuePointPatchField<vector>>(ppd);
            const pointField& pts = pvpd.patch().localPoints();
            pvpd.replace
            (
                0,
                pts.component(0)
              - sqrt(mag(1.0 - sqr(pts.component(1))))*sign(pts.component(0))
            );
            if (mesh.nGeometricD() > 2)
            {
                pvpd.replace
                (
                    2,
                    pts.component(2)
                  - sqrt(mag(1.0 - sqr(pts.component(1))))
                   *sign(pts.component(2))
                );
            }
        }
    }

    mesh.update();

    pointIOField& points0 =
        mesh.lookupObjectRef<pointIOField>("points0");

    runTime++;
    Info<<"Time = " << runTime.timeName() << nl <<endl;
    runTime.write();

    mesh.refine();
    forAll(pointDisplacement.boundaryField(), patchi)
    {
        pointPatchField<vector>& ppd =
            pointDisplacement.boundaryFieldRef()[patchi];
        if (isA<valuePointPatchField<vector>>(ppd))
        {
            valuePointPatchField<vector>& pvpd =
                dynamicCast<valuePointPatchField<vector>>(ppd);
            const vectorField pts
            (
                UIndirectList<vector>(points0, ppd.patch().meshPoints())()
            );
            pvpd.replace
            (
                0,
                pts.component(0)
              - sqrt(mag(1.0 - sqr(pts.component(1))))*sign(pts.component(0))
            );
            if (mesh.nGeometricD() > 2)
            {
                pvpd.replace
                (
                    2,
                    pts.component(2)
                  - sqrt(mag(1.0 - sqr(pts.component(1))))
                   *sign(pts.component(2))
                );
            }
        }
    }

    mesh.update();

    runTime++;
    Info<<"Time = " << runTime.timeName() << nl <<endl;
    runTime.write();

    mesh.refine();
    forAll(pointDisplacement.boundaryField(), patchi)
    {
        pointPatchField<vector>& ppd =
            pointDisplacement.boundaryFieldRef()[patchi];
        if (isA<valuePointPatchField<vector>>(ppd))
        {
            valuePointPatchField<vector>& pvpd =
                dynamicCast<valuePointPatchField<vector>>(ppd);
            const vectorField pts
            (
                UIndirectList<vector>(points0, ppd.patch().meshPoints())()
            );
            pvpd.replace
            (
                0,
                pts.component(0)
              - sqrt(mag(1.0 - sqr(pts.component(1))))*sign(pts.component(0))
            );
            if (mesh.nGeometricD() > 2)
            {
                pvpd.replace
                (
                    2,
                    pts.component(2)
                  - sqrt(mag(1.0 - sqr(pts.component(1))))
                   *sign(pts.component(2))
                );
            }
        }
    }

    mesh.update();

    runTime++;
    Info<<"Time = " << runTime.timeName() << nl <<endl;
    runTime.write();

    pointField points0Old(points0);

    mesh.movePoints(points0);

    runTime++;
    Info<<"Time = " << runTime.timeName() << nl <<endl;
    runTime.write();

    forAll(pointDisplacement.boundaryField(), patchi)
    {
        pointPatchField<vector>& ppd =
            pointDisplacement.boundaryFieldRef()[patchi];
        if (isA<valuePointPatchField<vector>>(ppd))
        {
            valuePointPatchField<vector>& pvpd =
                dynamicCast<valuePointPatchField<vector>>(ppd);
            const vectorField pts
            (
                UIndirectList<vector>(points0, ppd.patch().meshPoints())()
            );
            pvpd.replace
            (
                0,
                pts.component(0)
              - sqrt(mag(1.0 - sqr(pts.component(1))))
                *sign(pts.component(0))
            );
            if (mesh.nGeometricD() > 2)
            {
                pvpd.replace
                (
                    2,
                    pts.component(2)
                  - sqrt(mag(1.0 - sqr(pts.component(1))))
                   *sign(pts.component(2))
                );
            }
        }
    }
    pointDisplacement.boundaryFieldRef().updateCoeffs();
    mesh.update();

    runTime++;
    Info<<"Time = " << runTime.timeName() << nl <<endl;
    pointVectorField
    (
        "interp",
        volPointInterpolation::New(mesh).interpolate(cellDisplacement)
    ).write();
    runTime.write();

    e = 0;
    e[0] = 1;

    mesh.refine();

    runTime++;
    Info<<"Time = " << runTime.timeName() << nl <<endl;
    runTime.write();

    mesh.refine();

    runTime++;
    Info<<"Time = " << runTime.timeName() << nl <<endl;
    runTime.write();


    mesh.movePoints(points0);

    runTime++;
    Info<<"Time = " << runTime.timeName() << nl <<endl;
    runTime.write();


    Info<<"done"<<endl;

    return 0;
}
