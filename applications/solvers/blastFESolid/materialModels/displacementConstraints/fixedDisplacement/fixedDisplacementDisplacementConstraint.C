
#include "fixedDisplacementDisplacementConstraint.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace displacementConstraints
{
    defineTypeNameAndDebug(fixedDisplacement, 0);
    addToRunTimeSelectionTable
    (
        displacementConstraint,
        fixedDisplacement,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementConstraints::fixedDisplacement::fixedDisplacement
(
    const word& name,
    pointVectorField& D,
    pointVectorField& U,
    const dictionary& dict
)
:
    displacementConstraint(name, D, U, dict),
    displacement_
    (
        Function1<vector>::New
        (
            "displacement",
            D.time().userUnits(),
            dimLength,
            dict
        )
    )
{}


Foam::displacementConstraints::fixedDisplacement::fixedDisplacement
(
    const fixedDisplacement& fddc
)
:
    displacementConstraint(fddc),
    displacement_(fddc.displacement_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementConstraints::fixedDisplacement::~fixedDisplacement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::displacementConstraints::fixedDisplacement::constrain()
{
    const scalar t = D_.time().value();
    const scalar dt = D_.time().deltaTValue();
    const pointVectorField& D0 = D_.oldTime();

    const vector d = displacement_->value(t);
    forAll(nodes_, ni)
    {
        const label nodei = nodes_[ni];
        forAll(dims_, ci)
        {
            const label cmpti = dims_[ci];
            D_[nodei][cmpti] = d[cmpti];
            U_[nodei][cmpti] = (d[cmpti] - D0[nodei][cmpti])/dt;
        }
    }
}

// ************************************************************************* //
