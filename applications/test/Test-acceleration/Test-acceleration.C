/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fieldAccelerationSchemeList.H"
#include "mixedFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    runTime++;

    IOdictionary dict
    (
        IOobject
        (
            "regionProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ
        )
    );

    Info<< "Testing volScalarFields" << endl;
    fieldAccelerationSchemeList acceleration(mesh, dict);
    acceleration.addField(x, {0});
    acceleration.addField(y, {0});
    acceleration.addField(mixedx, {0});
    acceleration.addField(grady, {0});

    const label maxIter = dict.lookup<label>("maxIter");

    scalarField& px = x.boundaryFieldRef()[0];
    scalarField& prefValx = dynamicCast<mixedFvPatchScalarField>(mixedx.boundaryFieldRef()[0]).refValue();
    prefValx = px;
    scalarField& prefGradx =
        dynamicCast<mixedFvPatchScalarField>(mixedx.boundaryFieldRef()[0]).refGrad();
    prefGradx = 0.0;
    scalarField& pvalFracx =
        dynamicCast<mixedFvPatchScalarField>(mixedx.boundaryFieldRef()[0]).valueFraction();
    pvalFracx = 1.0;

    scalarField& py = y.boundaryFieldRef()[0];
    scalarField& pgrady =
        dynamicCast<fixedGradientFvPatchScalarField>(grady.boundaryFieldRef()[0]).gradient();
    pgrady = 0.0;

    mixedx.correctBoundaryConditions();
    grady.correctBoundaryConditions();

    for (label iter = 0; iter < maxIter; iter++)
    {
        acceleration.storePrevIter();

        px = pow3(py + 1.0);
        py = sqr(cos(px));
        prefGradx = py;
        prefValx = px;
        pvalFracx = 0.5;
        pgrady = 0.0;

        mixedx.correctBoundaryConditions();
        grady.correctBoundaryConditions();

        acceleration.relax(iter);

        Info<< "Iteration " << iter << endl;
        acceleration.print(Info);
        Info<< endl;

        if (acceleration.converged())
        {
            Info<<"converged" << nl << endl;
            break;
        }
    }
    Info<< "solution: " << 1.27574 << " " << 0.0845611 << endl;
    Info<< "computed: " << px[0] << " " << py[0] << nl << nl << endl;


    Info<< "Testing pointVectorFields" << endl;
    fieldAccelerationSchemeList pAcceleration(mesh, dict);
    pAcceleration.addField(pointx, {0});
    pAcceleration.addField(pointy, {0});

    vectorField& ppointx =
        dynamicCast<vectorField>(pointx.boundaryFieldRef()[0]);
    vectorField& ppointy =
        dynamicCast<vectorField>(pointy.boundaryFieldRef()[0]);
    for (label iter = 0; iter < maxIter; iter++)
    {
        pAcceleration.storePrevIter();

        forAll(ppointx, i)
        {
            ppointx[i] = cmptPow(ppointy[i] - 2.0*vector::one, vector::one*3);
            ppointy[i] = cmptPow(ppointx[i] + vector::one, vector::one*3);
        }

        pAcceleration.relax(iter);

        Info<< "Iteration " << iter << endl;
        pAcceleration.print(Info);
        Info<< endl;

        if (pAcceleration.converged())
        {
            Info<<"converged" << nl << endl;
            break;
        }
    }
    Info<< "solution: " << 0.398632 <<" "<< 2.73597 <<endl;
    Info<< "computed: " << ppointx[0][0] <<" "<<ppointy[0][0] << nl << endl;

    Info<< "Done" << endl;

    return 0;
}


// ************************************************************************* //
