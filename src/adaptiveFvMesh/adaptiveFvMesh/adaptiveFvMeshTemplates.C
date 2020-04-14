/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
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

#include "volMesh.H"
#include "fvPatchField.H"
#include "surfaceFields.H"

template <class T>
void Foam::adaptiveFvMesh::mapNewInternalFaces
(
    const labelList& faceMap
)
{
//     typedef GeometricField<T, fvsPatchField, surfaceMesh> GeoField;
//     HashTable<GeoField*> sFlds(this->objectRegistry::lookupClass<GeoField>());
//
//     const labelUList& owner = this->faceOwner();
//     const labelUList& neighbour = this->faceNeighbour();
//     const dimensionedScalar deltaN = 1e-8 / pow(average(this->V()), 1.0 / 3.0);
//
//     forAllIter(typename HashTable<GeoField*>, sFlds, iter)
//     {
//         GeoField& sFld = *iter();
//         if (mapSurfaceFields_.found(iter.key()))
//         {
//             if (debug)
//             {
//                 InfoInFunction << iter.key()<< endl;
//             }
//
//
//             // Create flat version of field
//             Field<T> tsFld(this->nFaces(), pTraits<T>::zero);
//             {
//                 forAll(sFld.internalField(), iFace)
//                 {
//                     tsFld[iFace] = sFld.internalField()[iFace];
//                 }
//                 forAll(sFld.boundaryField(), iPatch)
//                 {
//                     const fvsPatchField<T>& pfld = sFld.boundaryField()[iPatch];
//
//                     label start = pfld.patch().start();
//                     forAll(pfld, faceI)
//                     {
//                         tsFld[faceI+start] = pfld[faceI];
//                     }
//                 }
//             }
//
//
//             // Loop over all faces
//             for (label facei = 0; facei < nInternalFaces(); facei++)
//             {
//                 label oldFacei = faceMap[facei];
//
//                 //- map surface field on newly generated faces
//                 if (oldFacei == -1)
//                 {
//                     //- find owner and neighbour cell
//                     const cell& ownFaces = this->cells()[owner[facei]];
//                     const cell& neiFaces = this->cells()[neighbour[facei]];
//
//                     //- loop over all owner/neighbour cell faces
//                     //- and find already mapped ones (master-faces):
//                     T tmpValue = pTraits<T>::zero;
//                     scalar magFld = 0;
//                     label counter = 0;
//
//                     //- simple averaging of all neighbour master-faces
//                     forAll(ownFaces, iFace)
//                     {
//                         if (faceMap[ownFaces[iFace]] != -1)
//                         {
//                             tmpValue += tsFld[ownFaces[iFace]];
//                             magFld += mag(tsFld[ownFaces[iFace]]);
//                             counter++;
//                         }
//                     }
//
//                     forAll(neiFaces, iFace)
//                     {
//                         if (faceMap[neiFaces[iFace]] != -1)
//                         {
//                             tmpValue = tmpValue + tsFld[neiFaces[iFace]];
//                             magFld += mag(tsFld[neiFaces[iFace]]);
//                             counter++;
//                         }
//                     }
//
//                     if(counter > 0)
//                     {
//                         if
//                         (
//                             GeometricField<T, fvsPatchField, surfaceMesh>::typeName
//                                 == "surfaceScalarField"
//                         )
//                         {
//                             tmpValue /= counter;
//                         }
//                         else if
//                         (
//                             GeometricField<T, fvsPatchField, surfaceMesh>::typeName
//                                 == "surfaceVectorField"
//                         )
//                         {
//                             magFld /= counter;
//                             tmpValue *= magFld/(mag(tmpValue)+deltaN.value());
//                         }
//                         else
//                         {
//                             FatalErrorInFunction
//                                 << "mapping implementation only valid for"
//                                 << " scalar and vector fields! \n Field "
//                                 << sFld.name() << " is of type: "
//                                 << GeometricField<T, fvsPatchField, surfaceMesh>::typeName
//                                 << abort(FatalError);
//                         }
//                     }
//
//                     sFld[facei] = tmpValue;
//                 }
//             }
//         }
//     }
}


template<class Type>
void Foam::adaptiveFvMesh::correctBoundaries()
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    HashTable<GeoField*> flds(this->objectRegistry::lookupClass<GeoField>());

    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        GeoField& fld = *iter();

        //mimic "evaluate" but only for coupled patches (processor or cyclic)
        // and only for blocking or nonBlocking comms (no scheduled comms)
        if
        (
            Pstream::defaultCommsType == Pstream::commsTypes::blocking
         || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            label nReq = Pstream::nRequests();

            forAll(fld.boundaryField(), patchi)
            {
                if(fld.boundaryField()[patchi].coupled())
                {
                    fld.boundaryFieldRef()[patchi].initEvaluate
                    (
                        Pstream::defaultCommsType
                    );
                }
            }

            // Block for any outstanding requests
            if
            (
                Pstream::parRun()
             && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
            )
            {
                Pstream::waitRequests(nReq);
            }

            forAll(fld.boundaryField(), patchi)
            {
                if(fld.boundaryField()[patchi].coupled())
                {
                    fld.boundaryFieldRef()[patchi].evaluate
                    (
                        Pstream::defaultCommsType
                    );
                }
            }
        }
        else
        {
            //Scheduled patch updates not supported
            FatalErrorIn
            (
                "dynamicRefineBalancedFvMeshTemplates::correctBoundaries"
            )   << "Unsuported communications type "
                << Pstream::commsTypeNames[Pstream::defaultCommsType]
                << exit(FatalError);
        }
    }
}

