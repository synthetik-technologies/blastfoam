
#include "solidTractionBoundaryTraction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace boundaryTractions
{
    defineTypeNameAndDebug(solidTraction, 0);
    addToRunTimeSelectionTable
    (
        boundaryTraction,
        solidTraction,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryTractions::solidTraction::solidTraction
(
    const word& name,
    const dictionary& dict,
    const feMesh1& femesh,
    const pointVectorField* DPtr
)
:
    boundaryTraction(name, dict, femesh, DPtr),
    pressure_
    (
        Function1<scalar>::New
        (
            "pressure",
            femesh.time().userUnits(),
            dimPressure,
            dict
        )
    ),
    traction_
    (
        Function1<vector>::New
        (
            "traction",
            femesh.time().userUnits(),
            dimPressure,
            dict
        )
    )
{}


Foam::boundaryTractions::solidTraction::solidTraction
(
    const solidTraction& stbt
)
:
    boundaryTraction(stbt),
    pressure_(stbt.pressure_, false),
    traction_(stbt.traction_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boundaryTractions::solidTraction::~solidTraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::boundaryTractions::solidTraction::pressure
(
    const labelList& nodes,
    const scalarList& shape,
    const vector& n
) const
{
   return pressure_->value(mesh_.time().value());
}


Foam::vector Foam::boundaryTractions::solidTraction::traction
(
    const labelList& nodes,
    const scalarList& shape,
    const vector& n
) const
{
   return traction_->value(mesh_.time().value());
}


// ************************************************************************* //
