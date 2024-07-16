
#include "CONWEPPressureBoundaryTraction.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace boundaryTractions
{
    defineTypeNameAndDebug(CONWEPPressure, 0);
    addToRunTimeSelectionTable
    (
        boundaryTraction,
        CONWEPPressure,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryTractions::CONWEPPressure::CONWEPPressure
(
    const word& name,
    const dictionary& dict,
    const feMesh1& femesh,
    const pointVectorField* DPtr
)
:
    boundaryTraction(name, dict, femesh, DPtr),
    conwep_(dict),
    integrated_(dict.lookupOrDefault<bool>("integratedPressure", false)),
    pRef_(dict.lookup<scalar>("pRef"))
{}


Foam::boundaryTractions::CONWEPPressure::CONWEPPressure
(
    const CONWEPPressure& cpbt
)
:
    boundaryTraction(cpbt),
    conwep_(cpbt.conwep_),
    integrated_(cpbt.integrated_),
    pRef_(cpbt.pRef_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boundaryTractions::CONWEPPressure::~CONWEPPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::scalar Foam::boundaryTractions::CONWEPPressure::pressure
(
    const labelList& nodeLabels,
    const scalarList& shape,
    const vector& n
) const
{
    const UIndirectList<vector> nodes
    (
        femesh_.boundary()[patchID_].localNodes(),
        nodeLabels
    );

    const scalar t = mesh_.time().value();
    const scalar dt = mesh_.time().deltaTValue();
    const scalar t0 = t - dt;
    const vector nNeg(n);

    scalar p = 0.0;
    if (DPtr_)
    {
        const labelList& meshNodes =
            femesh_.boundary()[patchID_].meshNodes();
        const vectorField& D = *DPtr_;
        if (integrated_)
        {
            forAll(nodes, i)
            {
                p +=
                    shape[i]
                   *conwep_.impulse
                    (
                        nodes[i] + D[meshNodes[nodeLabels[i]]],
                        t0,
                        t,
                        nNeg
                    )/dt;
            }
        }
        else
        {
            forAll(nodes, i)
            {
                p +=
                    shape[i]
                   *conwep_.overpressure
                    (
                        nodes[i] + D[meshNodes[nodeLabels[i]]],
                        t,
                        nNeg
                    );
            }
        }
    }
    else
    {
        if (integrated_)
        {
            forAll(nodes, i)
            {
                p +=
                    shape[i]
                   *conwep_.impulse(nodes[i], t0, t, nNeg)/dt;
            }
        }
        else
        {
            forAll(nodes, i)
            {
                p +=
                    shape[i]
                   *conwep_.overpressure(nodes[i], t, nNeg);
            }
        }
    }
    return p + pRef_;
}


Foam::vector Foam::boundaryTractions::CONWEPPressure::traction
(
    const labelList& nodeLabels,
    const scalarList& shape,
    const vector& n
) const
{
   return Zero;
}


// ************************************************************************* //
