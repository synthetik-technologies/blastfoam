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

#include "boundingBoxLimiter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(boundingBoxLimiter, 0);
    addToRunTimeSelectionTable
    (
        crackPathLimiter, boundingBoxLimiter, dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::boundingBoxLimiter::calcFacesAllowedToBreak() const
{
    if (facesAllowedToBreakPtr_)
    {
        FatalErrorIn
        (
            "void Foam::boundingBoxLimiter::calcFacesAllowedToBreak() const"
        ) << "pointer already set" << abort(FatalError);
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

    Info<< nl << "There are " << gSum(facesAllowedToBreak.internalField())
        << " potential internal crack faces" << nl
        << "There are " << gSum(facesAllowedToBreak.boundaryField())/2
        << " potential coupled boundary crack faces" << endl;


    // It is currently not possible to directly visualise surface fields in
    // ParaView, so we create a volume field to show cells adjacent to potential
    // cohesive faces

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

    forAll(facesAllowedToBreakI, faceI)
    {
        if (facesAllowedToBreakI[faceI] > SMALL)
        {
            crackLimiterBoxesI[owner[faceI]] = 1.0;
            crackLimiterBoxesI[neighbour[faceI]] = 1.0;
        }
    }

    forAll(crackLimiterBoxes.boundaryField(), patchI)
    {
        crackLimiterBoxes.boundaryFieldRef()[patchI] =
            facesAllowedToBreak.boundaryField()[patchI];
    }

    if (mesh.time().outputTime())
    {
        Info<< "Writing cohesiveZone field" << endl;
        crackLimiterBoxes.write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from dictionary
Foam::boundingBoxLimiter::boundingBoxLimiter
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


Foam::boundingBoxLimiter::~boundingBoxLimiter()
{
    clearOut();
}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

const Foam::surfaceScalarField&
Foam::boundingBoxLimiter::facesAllowedToBreak() const
{
    if (!facesAllowedToBreakPtr_)
    {
        calcFacesAllowedToBreak();
    }

    return *facesAllowedToBreakPtr_;
}


void Foam::boundingBoxLimiter::clearOut()
{
    deleteDemandDrivenData(facesAllowedToBreakPtr_);
}


// ************************************************************************* //
