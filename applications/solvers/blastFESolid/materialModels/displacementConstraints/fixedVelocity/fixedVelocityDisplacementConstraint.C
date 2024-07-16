
#include "fixedVelocityDisplacementConstraint.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace displacementConstraints
{
    defineTypeNameAndDebug(fixedVelocity, 0);
    addToRunTimeSelectionTable
    (
        displacementConstraint,
        fixedVelocity,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementConstraints::fixedVelocity::fixedVelocity
(
    const word& name,
    pointVectorField& D,
    pointVectorField& U,
    const dictionary& dict
)
:
    displacementConstraint(name, D, U, dict),
    velocity_(Function1<vector>::New("velocity", dict))
{}


Foam::displacementConstraints::fixedVelocity::fixedVelocity
(
    const fixedVelocity& fvdc
)
:
    displacementConstraint(fvdc),
    velocity_(fvdc.velocity_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementConstraints::fixedVelocity::~fixedVelocity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::displacementConstraints::fixedVelocity::constrain()
{
    const scalar t = D_.time().value();
    const scalar dt = D_.time().deltaTValue();
    const pointVectorField& D0 = D_.oldTime();

    const vector u = velocity_->value(t);
    forAll(nodes_, ni)
    {
        const label nodei = nodes_[ni];
        forAll(dims_, ci)
        {
            const label cmpti = dims_[ci];
            D_[nodei][cmpti] = D0[nodei][cmpti] + dt*u[cmpti];
            U_[nodei][cmpti] = u[cmpti];
        }
    }
}

// ************************************************************************* //
