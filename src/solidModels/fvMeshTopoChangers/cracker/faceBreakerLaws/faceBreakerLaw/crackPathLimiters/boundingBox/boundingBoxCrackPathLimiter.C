/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "boundingBoxCrackPathLimiter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace crackPathLimiters
{
    defineTypeNameAndDebug(boundingBox, 0);
    addToRunTimeSelectionTable(crackPathLimiter, boundingBox, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::crackPathLimiters::boundingBox::calcFacesAllowedToBreak() const
{
    if (facesAllowedToBreakPtr_)
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    const fvMesh& mesh = this->mesh();

    facesAllowedToBreakPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "facesAllowedToBreak",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        );

    surfaceScalarField& facesAllowedToBreak = *facesAllowedToBreakPtr_;

    const vectorField& CfI = mesh.Cf().internalField();
    scalarField& facesAllowedToBreakI =
        facesAllowedToBreak.primitiveFieldRef();

    forAll(boundingBoxes_, boxI)
    {
        forAll(facesAllowedToBreakI, faceI)
        {
            if (boundingBoxes_[boxI].contains(CfI[faceI]))
            {
                facesAllowedToBreakI[faceI] = 1.0;
            }
        }

        forAll(facesAllowedToBreak.boundaryField(), patchI)
        {
            if (mesh.boundaryMesh()[patchI].coupled())
            {
                scalarField& pFacesAllowedToBreak =
                    facesAllowedToBreak.boundaryFieldRef()[patchI];
                const vectorField& pCf =
                    mesh.Cf().boundaryField()[patchI];

                forAll(pFacesAllowedToBreak, faceI)
                {
                    if (boundingBoxes_[boxI].contains(pCf[faceI]))
                    {
                        pFacesAllowedToBreak[faceI] = 1.0;
                    }
                }
            }
        }
    }

    DebugInfo
        << nl
        << "There are " << gSum(facesAllowedToBreak.internalField())
        << " potential internal crack faces" << nl
        << "There are " << gSum(facesAllowedToBreak.boundaryField())/2
        << " potential coupled boundary crack faces" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from dictionary
Foam::crackPathLimiters::boundingBox::boundingBox
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    crackPathLimiter(name, mesh, dict),
    facesAllowedToBreakPtr_(NULL),
    boundingBoxes_(dict.lookup("boundingBoxes"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::crackPathLimiters::boundingBox::~boundingBox()
{
    clearOut();
}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::crackPathLimiters::boundingBox::facesAllowedToBreak() const
{
    if (!facesAllowedToBreakPtr_)
    {
        calcFacesAllowedToBreak();
    }

    return *facesAllowedToBreakPtr_;
}


void Foam::crackPathLimiters::boundingBox::clearOut()
{
    deleteDemandDrivenData(facesAllowedToBreakPtr_);
}


bool Foam::crackPathLimiters::boundingBox::write() const
{
    // It is currently not possible to directly visualise surface fields in
    // ParaView, so we create a volume field to show cells adjacent to potential
    // cohesive faces
    const fvMesh& mesh = this->mesh();
    volScalarField crackLimiterBoxes
    (
        IOobject
        (
            "crackLimiterBoxes",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    scalarField& crackLimiterBoxesI =
        crackLimiterBoxes.primitiveFieldRef();
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const surfaceScalarField& facesAllowedToBreak =
        this->facesAllowedToBreak();
    forAll(facesAllowedToBreak, faceI)
    {
        if (facesAllowedToBreak[faceI] > SMALL)
        {
            crackLimiterBoxesI[owner[faceI]] = 1.0;
            crackLimiterBoxesI[neighbour[faceI]] = 1.0;
        }
    }

    volScalarField::Boundary& bcrackLimiterBoxes =
        crackLimiterBoxes.boundaryFieldRef();
    const surfaceScalarField::Boundary& bfacesAllowedToBreak =
        facesAllowedToBreak.boundaryField();
    forAll(bcrackLimiterBoxes, patchI)
    {
        bcrackLimiterBoxes[patchI] =
            bfacesAllowedToBreak[patchI];
    }

    DebugInfo<< "Writing cohesiveZone field" << endl;
    return crackLimiterBoxes.write();
}

// ************************************************************************* //
