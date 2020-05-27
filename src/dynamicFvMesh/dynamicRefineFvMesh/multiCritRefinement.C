/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "multiCritRefinement.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::multiCritRefinement::readMultiCritRefinementDict()
{
    dictionary dynamicMeshDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        )
    );

    nBufferLayers_ = readLabel(dynamicMeshDict.subDict("dynamicRefineFvMeshCoeffs").lookup("nBufferLayers"));

    if( dynamicMeshDict.isDict("multiCritRefinementControls") )
    {
        dictionary refineControlDict =
            dynamicMeshDict.subDict("multiCritRefinementControls");

        enableMultiCritRefinementControl_ =
            Switch(refineControlDict.lookup("enableMultiCritRefinementControl"));

        if( enableMultiCritRefinementControl_ )
        {

            dictionary multiCritRefinementControls =
                dynamicMeshDict.subDict("multiCritRefinementControls");


            // Overwrite field name entry in dynamicRefineFvMeshCoeffs?
            // For now you just have to be smart and enter
            // 'multiCritRefinementField' for the field name manually

            // Read HashTable of field-refinement scalars
            if( multiCritRefinementControls.found("fields") )
            {
                fields_ = HashTable< dictionary >
                (
                    multiCritRefinementControls.lookup("fields")
                );
            }

            // Read HashTable of gradient-refinement scalars
            if( multiCritRefinementControls.found("gradients") )
            {
                gradFields_ = HashTable< dictionary >
                (
                    multiCritRefinementControls.lookup("gradients")
                );
            }

            // Read HashTable of curl-refinement vectors
            if( multiCritRefinementControls.found("curls") )
            {
                curlFields_ = HashTable< dictionary >
                (
                    multiCritRefinementControls.lookup("curls")
                );
            }

            // Read interface refinement data
            if( multiCritRefinementControls.found("interface") )
            {
                interface_ = HashTable< dictionary >
                (
                    multiCritRefinementControls.lookup("interface")
                );
            }

            // Read refinement regions
            if( multiCritRefinementControls.found("regions") )
            {
                refinedRegions_ = PtrList<entry>
                (
                    multiCritRefinementControls.lookup("regions")
                );
            }

            return true;
        }else{
            return false;
        }
    }
    else
    {
        return false;
    }
}


void Foam::multiCritRefinement::applyCritEntries(word critType, dictionary critDict, word critName)
{

    scalar minValue = readScalar(critDict.lookup("minValue"));
    scalar maxValue = readScalar(critDict.lookup("maxValue"));
    scalar refineLevel = readScalar(critDict.lookup("refineLevel"));

    //- allow different criteria applied on one field, using different hash keys.
    word fldName= "";
    if (critDict.found("fieldName"))
    {
        fldName = word(critDict.lookup("fieldName"));
    } else {
        fldName = critName;
    }

    Field<scalar> refFld(mesh_.nCells(), 0.0);

    if (critType == "field")
    {
        //- get a handle of the field
        const volScalarField& fld = mesh_.lookupObject<volScalarField>(fldName);

        refFld = fld;
    }

    if (critType == "grad")
    {
        //- get a handle of the field
        const volVectorField& fld = mesh_.lookupObject<volVectorField>(fldName);

        Field<scalar> cubeRtV = Foam::pow(mesh_.V(),1.0/3.0);

        refFld = mag(fvc::grad(fld)) * cubeRtV;
    }

    if (critType == "curl")
    {
        //- get a handle of the field
        const volVectorField& fld = mesh_.lookupObject<volVectorField>(fldName);

        refFld = mag(fvc::curl(fld));
    }

    // Limit the value of refFld based on its max level
    forAll(refFld, cellI)
    {

        if ((refFld[cellI] >= minValue) && (refFld[cellI] <= maxValue))
        {
            // increase targetLevel up to refineLevel
            // BUT: do not decrease if cell already marked for higher refinement level by previous criterion
            targetLevel_[cellI] = max(targetLevel_[cellI], refineLevel);
        }
    }

    //- AddLayer Keyword
    scalar nAddLayers(0);
    if (critDict.found("nAddLayers"))
    {
        nAddLayers = readScalar(critDict.lookup("nAddLayers"));
    }

    if (nAddLayers > 0)
    {
        //- create a volField of the labelList targetLevel
            // Set the internal refinement field to zero to start with
        volScalarField& tmp = *multiCritRefinementFieldPtr_;
        tmp = dimensionedScalar("zero",dimless,0.0);
        volScalarField tLevel(tmp * 0.0);
        forAll(targetLevel_, cellI)
        {
            tLevel[cellI] = targetLevel_[cellI];
        }

        // #nAddLayers cells per additional level
        for(label j=0; j < nAddLayers; j++)
        {
            //- select the area with targetLevel==refineLevel
            volScalarField finest(pos(tLevel - refineLevel + SMALL));

            //- add +1 to targetLevel on the enlarged stencil
            tLevel += pos( fvc::average(fvc::interpolate(finest) - SMALL)) - finest;
        }

        //- copy the new results onto the target labelList
        forAll(targetLevel_, cellI)
        {
            targetLevel_[cellI] = tLevel[cellI];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiCritRefinement::multiCritRefinement(
    const labelList& cellLevel,
    const fvMesh& mesh
)
:
    cellLevel_(cellLevel),
    mesh_(mesh),
    multiCritRefinementFieldPtr_(NULL),
    targetLevelPtr_(NULL),
    targetLevel_(mesh.nCells(),0),
    isLevelPtr_(NULL),
    fields_(),
    gradFields_(),
    curlFields_(),
    interface_(),
    refinedRegions_(),
    enableMultiCritRefinementControl_(false),
    nBufferLayers_(0)
{
    // Read static part of dictionary
    readMultiCritRefinementDict();

    multiCritRefinementFieldPtr_ = new volScalarField
        (
            IOobject
            (
                "multiCritRefinementField",
                mesh.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        );

        targetLevelPtr_ = new volScalarField
        (
            IOobject
            (
                "targetLevel",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        );

        isLevelPtr_ = new volScalarField
        (
            IOobject
            (
                "isLevel",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        );

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiCritRefinement::~multiCritRefinement()
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiCritRefinement::updateRefinementField()
{

    if( !readMultiCritRefinementDict() )
    {
        return;
    }

    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict("dynamicRefineFvMeshCoeffs")
    );
    scalar globalMaxRefLevel = readScalar(refineDict.lookup("maxRefinement"));

    Info<< "Calculating multiCritRefinement field" << endl;

    volScalarField& intRefFld = *multiCritRefinementFieldPtr_;
    volScalarField& targetFld = *targetLevelPtr_;
    volScalarField& currFld = *isLevelPtr_;

    // Set the internal refinement field to zero to start with
    intRefFld = dimensionedScalar("zero",dimless,0.0);

    // init list for refine/unrefine field (-1=unrefine, 1=refine, 0=do nothing)
    labelList markRefineCells (mesh_.nCells(), 0);

    // init list for target refinement level per cell
    targetLevel_ = markRefineCells;

    Field<scalar> refFld(mesh_.nCells(),0.0);

    // First fields
    {
        List<word> fieldNames = fields_.toc();


        forAll(fieldNames, i)
        {

            word fldEntry = fieldNames[i];

            //- read criteria dict entries and calculate target level
            applyCritEntries("field", fields_[fldEntry], fldEntry);
        }
    }

    // Then gradients
    {
        List<word> gradFieldNames = gradFields_.toc();

        forAll(gradFieldNames, i)
        {
             word fldEntry = gradFieldNames[i];

            //- read criteria dict entries and calculate target level
            applyCritEntries("grad", gradFields_[fldEntry], fldEntry);
        }
    }

    // Then curls
    {
        List<word> curlFieldNames = curlFields_.toc();

        forAll(curlFieldNames, i)
        {
            word fldEntry = curlFieldNames[i];

            //- read criteria dict entries and calculate target level
            applyCritEntries("curl", curlFields_[fldEntry], fldEntry);
        }
    }

    // Then at the interface (assumed to be always the maximum refinement level)
    {
        List<word> interfaceRefineField = interface_.toc();
        forAll(interfaceRefineField, i)
        {
            word fldName = interfaceRefineField[i];

            // read region of maximum refinement levels inside and outside of interface indicator field
            // (does not need to be alpha, can also be a concentration field)
            scalar innerRefLayers = readScalar(interface_[fldName].lookup("innerRefLayers"));
            scalar outerRefLayers = readScalar(interface_[fldName].lookup("outerRefLayers"));

            scalar nAddLayers(0);
            if (interface_[fldName].found("nAddLayers"))
            {
                nAddLayers = readScalar(interface_[fldName].lookup("nAddLayers"));
            }
            label refineLevel = globalMaxRefLevel;

            if (interface_[fldName].found("maxRefineLevel"))
            {
                refineLevel = readScalar(interface_[fldName].lookup("maxRefineLevel"));

                // to avoid user input mistakes, limit the value with the maximum allowed
                refineLevel = min(globalMaxRefLevel, refineLevel);
            }

            const volScalarField& fld = mesh_.lookupObject<volScalarField>(fldName);

            volScalarField isInterface(intRefFld * 0.0);
            isInterface = dimensionedScalar("isInterface",dimless,0.0);

            surfaceScalarField deltaAlpha(mag(fvc::snGrad(fld) / mesh_.deltaCoeffs()));

            const unallocLabelList& owner = mesh_.owner();
            const unallocLabelList& neighbour = mesh_.neighbour();

            forAll(deltaAlpha, faceI)
            {
            //TODO: add a min max value of the specified field to be marked for refinement
            // NOW: hard-coded method snGradAlpha: checks mag(alpha_N - alpha_P) > 0.1 ? 1:0
                if (deltaAlpha[faceI] > 0.1) // currently a fixed prescribed value; should be read from díctionary
                {
                    label own = owner[faceI];
                    label nei = neighbour[faceI];

                    // set isInterface field to one
                    isInterface[own] = 1.0;
                    isInterface[nei] = 1.0;
                }
            }

            // assumed fld=0.5*(fldMax+fldMin) defines the interface
            dimensionedScalar fldInterfaceValue(0.5*(gMax(fld)+gMin(fld)));

            // Face-wise interplation for 2D meshes since point interpolation does not work
            // on empty patches
            if (mesh_.nGeometricD() == 2)
            {

                //-DD: implementation based on face interpolation
                //     which results in slower transport in diagonal direction
                // add inner refinement layers
                for(label i=0; i < innerRefLayers; i++)
                {
                   isInterface += neg(- fvc::average(fvc::interpolate(isInterface)) * pos(fld - fldInterfaceValue));
                   isInterface = neg(- isInterface);
                }

                // add outer refinement layers
                for(label i=0; i < outerRefLayers; i++)
                {
                   isInterface += neg(- fvc::average(fvc::interpolate(isInterface)) * pos(fldInterfaceValue - fld));
                   isInterface = neg(- isInterface);
                }

                forAll(isInterface, cellI)
                {
                   if (isInterface[cellI] > 0.5)
                   {
                       targetLevel_[cellI] = max(targetLevel_[cellI], refineLevel);
                   }
                }
            } else {

                //-DD: version using volPointInterpolation (direction independent buffer layer)
                const volPointInterpolation& pInterp = volPointInterpolation::New(mesh_);
                //const fvMesh& mesh = fld.mesh();

                // add inner refinement layers
                for(label i=0; i < innerRefLayers; i++)
                {
                    volScalarField markInner(isInterface*pos(fld - fldInterfaceValue));
                    pointScalarField markLayerP(pInterp.interpolate(markInner));

                    forAll(mesh_.C(), cellI)
                    {
                        scalar sum = 0.;
                        label nPoints = 0;

                        forAll(mesh_.cellPoints()[cellI], pointI)
                        {
                            sum += markLayerP[mesh_.cellPoints()[cellI][pointI]];
                            nPoints++;
                        }
                        if (nPoints > 0)
                        {
                            sum /= nPoints;
                        }
                        isInterface[cellI] += sum;
                    }
                }
                isInterface = pos(isInterface - SMALL);

                // add outer refinement layers
                for(label i=0; i < outerRefLayers; i++)
                {
                    volScalarField markOuter(isInterface*pos(fldInterfaceValue - fld));
                    pointScalarField markLayerP(pInterp.interpolate(markOuter));

                    forAll(mesh_.C(), cellI)
                    {
                        scalar sum = 0.;
                        label nPoints = 0;

                        forAll(mesh_.cellPoints()[cellI], pointI)
                        {
                            sum += markLayerP[mesh_.cellPoints()[cellI][pointI]];
                            nPoints++;
                        }
                        if (nPoints > 0)
                        {
                            sum /= nPoints;
                        }
                        isInterface[cellI] += sum;
                    }
                }
                isInterface = pos(isInterface - SMALL);

                forAll(isInterface, cellI)
                {
                    if (isInterface[cellI] > 0.5)
                    {
                        targetLevel_[cellI] = max(targetLevel_[cellI], refineLevel);
                    }
                }
            }

            //-DD: old implementation based on face interpolation
            //     which results in slower transport in diagonal direction
            // expand additional layers if specified:
            if (nAddLayers > 0)
            {
               // loop over additional layers
               for(label i=1; i < refineLevel; i++)
               {
                   // #nAddLayers cells per additional level
                   for(label j=0; j < nAddLayers; j++)
                   {
                       isInterface += neg(- fvc::average(fvc::interpolate(isInterface)));
                       isInterface = neg(- isInterface);
                   }

                   forAll(isInterface, cellI)
                   {
                       if (isInterface[cellI] == 1.0)
                       {
                           targetLevel_[cellI] = max(targetLevel_[cellI], (refineLevel - i));
                       }
                   }
               }
            }

        }
    }

    //- Finally geometric regions (force the mesh to stay refined near key features)
    {
        forAll(refinedRegions_, regionI)
        {
            const entry& region = refinedRegions_[regionI];

            autoPtr<topoSetSource> source =
                topoSetSource::New(region.keyword(), mesh_, region.dict());

            cellSet selectedCellSet
            (
                mesh_,
                "cellSet",
                mesh_.nCells()/10+1 //estimate
            );

            source->applyToSet
            (
                topoSetSource::NEW,
                selectedCellSet
            );

            const labelList cells = selectedCellSet.toc();

            label minLevel = readLabel(region.dict().lookup("minLevel"));

            forAll(cells, i)
            {
                const label& cellI = cells[i];

                targetLevel_[cellI] = max(targetLevel_[cellI], minLevel);
            }
        }
    }

    //-DD: buffer layer based on targetLevel field to prevent 2-to-1 refinement
    {
        volScalarField blockedLevel(targetFld*0.0);

        for (label currLayer=globalMaxRefLevel; currLayer>=1; currLayer--)
        {
            forAll (targetLevel_, cellI)
            {
                if (targetLevel_[cellI] >= currLayer)
                {
                    blockedLevel[cellI] = targetLevel_[cellI];
                }
            }

            for (label i=0; i<nBufferLayers_; i++)
            {
                blockedLevel = max(blockedLevel, pos(fvc::average(fvc::interpolate(blockedLevel)) - SMALL)*(currLayer-1));
            }

            labelList blockLev(blockedLevel.internalField().size(), 0);
            forAll (blockLev, cellI)
            {
                blockLev[cellI] = blockedLevel[cellI];
            }
            targetLevel_ = max(targetLevel_, blockLev);
        }
    }

    // mark cells to be refined/unrefined based on above criteria:
    // simply check, if targetLevel lower or higher cellLevel
    forAll (intRefFld.internalField(), cellI)
    {
        intRefFld.primitiveFieldRef()[cellI] = targetLevel_[cellI] - cellLevel_[cellI];
        targetFld.primitiveFieldRef()[cellI] = targetLevel_[cellI];
        currFld.primitiveFieldRef()[cellI] =  cellLevel_[cellI];
    }

    intRefFld.correctBoundaryConditions();

    Info<<"Min,max refinement field = " << Foam::min(intRefFld).value() << ", "
        << Foam::max(intRefFld).value() << endl;
}


// ************************************************************************* //
