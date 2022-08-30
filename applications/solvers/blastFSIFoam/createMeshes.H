regionProperties rp(runTime);
IOdictionary regionPropertiesDict
(
    IOobject
    (
        "regionProperties",
        runTime.time().constant(),
        runTime.db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

#include "createFluidMeshes.H"
#include "createSolidMeshes.H"
