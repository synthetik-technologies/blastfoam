#include "readGravitationalAcceleration.H"

Info<< "Creating twoPhaseSystem\n" << endl;
phaseSystem fluid(mesh);
autoPtr<timeIntegrator> integrator(timeIntegrator::New(mesh));
integrator->addSystem(fluid);

PtrList<phaseModel>& phases(fluid.phases());

autoPtr<errorEstimator> error;
if (!isA<staticFvMesh>(mesh))
{
    error =
        errorEstimator::New
        (
            mesh,
            IOdictionary
            (
                IOobject
                (
                    "dynamicMeshDict",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
}
