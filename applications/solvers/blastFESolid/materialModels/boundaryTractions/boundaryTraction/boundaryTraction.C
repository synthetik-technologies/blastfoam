
#include "boundaryTraction.H"
#include "faceSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(boundaryTraction, 0);
    defineRunTimeSelectionTable(boundaryTraction, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryTraction::boundaryTraction
(
    const word& name,
    const dictionary& dict,
    const feMesh1& femesh,
    const pointVectorField* DPtr
)
:
    name_(name),
    femesh_(femesh),
    mesh_(femesh.mesh()),
    DPtr_(DPtr),
    patchID_(mesh_.boundaryMesh()[name].index())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boundaryTraction::~boundaryTraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::boundaryTraction::getTraction
(
    const labelList& nodeLabels,
    const scalarList& shape,
    const vector& n
) const
{
    return
        this->traction(nodeLabels, shape, n)
      - n*this->pressure(nodeLabels, shape, n);
}
// ************************************************************************* //

