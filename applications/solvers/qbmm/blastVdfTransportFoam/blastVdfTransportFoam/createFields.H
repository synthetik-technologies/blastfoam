// Read gravitational acceleration
Info<< "\nReading g" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedVector(dimAcceleration, Zero)
);

Info<< "Reading phaseProperties dictionary\n" << endl;
IOdictionary phaseProperties
(
     IOobject
     (
          "phaseProperties",
          runTime.constant(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
     )
);

coupledMultiphaseCompressibleSystem fluid(mesh, phaseProperties);
autoPtr<timeIntegrator> integrator(timeIntegrator::New(mesh));
integrator->addSystem(fluid);

volScalarField& p = fluid.p();
const volVectorField& U = fluid.U();
const volScalarField& rho = fluid.rho();
const surfaceScalarField& phi = fluid.phi();
fluid.update();


Info<< "Reading populationBalanceProperties\n" << endl;

IOdictionary populationBalanceProperties
(
    IOobject
    (
        "populationBalanceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
word phaseName(populationBalanceProperties.lookup("phaseName"));
dimensionedScalar residualAlpha
(
    "residualAlpha",
    dimless,
    populationBalanceProperties
);
dimensionedScalar rhop
(
    "rho",
    dimDensity,
    populationBalanceProperties
);
word dragModel(populationBalanceProperties.lookup("dragModel"));

volVectorField Up
(
    IOobject
    (
        IOobject::groupName("U", phaseName),
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimVelocity, Zero)
);

autoPtr<ODEPopulationBalanceModel> populationBalance
(
    ODEPopulationBalanceModel::New
    (
        phaseName, populationBalanceProperties, phi
    )
);
integrator->addSystem(populationBalance());

velocityQuadratureApproximation& quadrature =
    mesh.lookupObjectRef<velocityQuadratureApproximation>
    (
        IOobject::groupName
        (
            "quadratureProperties",
            phaseName
        )
    );

const labelListList& momentOrders = quadrature.momentOrders();
volVelocityMomentFieldSet& moments = quadrature.moments();
mappedPtrList<volVelocityNode>& nodes = quadrature.nodes();
const labelList& velocityIndexes = quadrature.nodes()[0].velocityIndexes();
const label sizeIndex = nodes[0].sizeIndex();
const label nDims = momentOrders[0].size();

bool computeVariance = false;
forAll(momentOrders, mi)
{
    forAll(velocityIndexes, cmpt)
    {
        if (momentOrders[mi][velocityIndexes[cmpt]] >= 2)
        {
            computeVariance = true;
        }
    }
}

const volScalarField& alpha = moments(labelList(nDims, 0));

autoPtr<volTensorField> Sigmap;
autoPtr<volScalarField> Thetap;

if (computeVariance)
{
    Sigmap.set
    (
        new volTensorField
        (
            IOobject
            (
                IOobject::groupName("Sigma", phaseName),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedTensor("zero", sqr(dimVelocity), Zero)
        )
    );

    Thetap.set
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("Theta", phaseName),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", sqr(dimVelocity), 0.0)
        )
    );
}

fluid.setDispersedVolumeFraction(alpha);
fluid.decode();
fluid.encode();
#include "computeParticleFields.H"
