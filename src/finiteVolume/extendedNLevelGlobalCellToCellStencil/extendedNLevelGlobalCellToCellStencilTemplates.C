/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "extendedNLevelGlobalCellToCellStencil.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class StencilType>
template<class Type, template<class> class ListType>
void Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::collectData
(
    const UList<Type>& fld,
    List<ListType<Type>>& stencilFld
) const
{
    if (!mapPtr_.valid())
    {
        update();
    }

    // 1. Construct cell data in compact addressing
    List<Type> flatFld(map().constructSize(), Zero);

    // Insert my internal values
    forAll(fld, celli)
    {
        flatFld[celli] = fld[celli];
    }

    // Do all swapping
    map().distribute(flatFld);

    // 2. Pull to stencil
    stencilFld.setSize(cellCells_.size());

    forAll(cellCells_, celli)
    {
        const labelList& compactCells = cellCells_[celli];

        stencilFld[celli].setSize(compactCells.size());

        forAll(compactCells, i)
        {
            stencilFld[celli][i] = flatFld[compactCells[i]];
        }
    }
}


template<class StencilType>
template<class Type>
void
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::collectOwnerData
(
    const UList<Type>& fld,
    Map<Type>& mapFld
) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    DynamicList<label> requests(nonlocalOwners_.size());
    DynamicList<label> requestedCells(nonlocalOwners_.size());
    forAllConstIter(Map<label>, nonlocalOwners_, iter)
    {
        requests.append(iter());
        requestedCells.append(iter.key());
    }

    // Needed for reverseDistribute
    label constructSize = requests.size();

    // make the map
    autoPtr<mapDistribute> map(buildMap(requests));

    // Send requests
    labelList sendCells(requestedCells);
    map().distribute(sendCells);

    // Add neighbours to be sent using the ordering provided
    List<Type> sendData(sendCells.size());
    forAll(sendCells, i)
    {
        sendData[i] = fld[gIndexPtr_->toLocal(sendCells[i])];
    }

    // Send the data back
    map().reverseDistribute(constructSize, sendData);

    forAll(requestedCells, i)
    {
        mapFld.insert(requestedCells[i], sendData[i]);
    }
}


template<class StencilType>
template<class Type>
void
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::collectOwnerData
(
    Map<Type>& mapFld
) const
{
    DynamicList<label> requests(nonlocalOwners_.size());
    DynamicList<label> requestedCells(nonlocalOwners_.size());
    forAllConstIter(Map<label>, nonlocalOwners_, iter)
    {
        requests.append(iter());
        requestedCells.append(iter.key());
    }

    // Needed for reverseDistribute
    label constructSize = requests.size();

    // make the map
    autoPtr<mapDistribute> map(buildMap(requests));

    // Send requests
    labelList sendCells(requestedCells);
    map().distribute(sendCells);

    // Add neighbours to be sent using the ordering provided
    List<Type> sendData(sendCells.size());
    forAll(sendCells, i)
    {
        sendData[i] = mapFld[sendCells[i]];
    }

    // Send the data back
    map().reverseDistribute(constructSize, sendData);

    forAll(requestedCells, i)
    {
        mapFld.insert(requestedCells[i], sendData[i]);
    }
}


template<class StencilType>
template<class Type>
void
Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::collectNbrData
(
    Map<Type>& mapFld
) const
{
    DynamicList<label> requests(nonlocalCells_.size());
    DynamicList<label> requestedCells(nonlocalCells_.size());
    forAllConstIter(Map<label>, nonlocalCells_, iter)
    {
        requests.append(iter());
        requestedCells.append(iter.key());
    }

    // Needed for reverseDistribute
    label constructSize = requests.size();

    // make the map
    autoPtr<mapDistribute> map(buildMap(requests));

    // Send requests
    labelList sendCells(requestedCells);
    map().distribute(sendCells);

    // Add neighbours to be sent using the ordering provided
    List<Type> sendData(sendCells.size());
    forAll(sendCells, i)
    {
        sendData[i] = mapFld[sendCells[i]];
    }

    // Send the data back
    map().reverseDistribute(constructSize, sendData);

    forAll(requestedCells, i)
    {
        mapFld.insert(requestedCells[i], sendData[i]);
    }
}


template<class StencilType>
template<class Type, class BinaryOp>
void Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::reduce
(
    const Map<Type>& mapFld,
    UList<Type>& fld,
    const BinaryOp& bop
) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    // Create the requests for new cell neighbours
    DynamicList<label> requests(fld.size());
    DynamicList<label> sendCells(fld.size());
    DynamicList<Type> sendData(fld.size());
    forAllConstIter
    (
        Map<cellStencil>,
        cellCellMap(),
        iter
    )
    {
        if (!iter().isLocal())
        {
            const label proci = gIndexPtr_->whichProcID(iter().owner());
            requests.append(proci);
            sendCells.append(gIndexPtr_->toLocal(proci, iter().owner()));
            sendData.append(mapFld[iter().owner()]);
        }
    }

    // make the map
    autoPtr<mapDistribute> map(buildMap(requests));

    // Send requests
    map().distribute(sendCells);
    map().distribute(sendData);

    // Add neighbours to be sent using the ordering provided
    forAll(sendCells, i)
    {
        const label celli = sendCells[i];
        fld[celli] = bop(fld[celli], sendData[i]);
    }
}


template<class StencilType>
template<class Type, class WeightType, template<class> class ListType>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<WeightType, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
> Foam::extendedNLevelGlobalCellToCellStencil<StencilType>::weightedSum
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const UList<ListType<WeightType>>& stencilWeights
) const
{
    typedef typename outerProduct<WeightType, Type>::type WeightedType;
    typedef GeometricField<WeightedType, fvPatchField, volMesh>
        WeightedFieldType;

    const fvMesh& mesh = fld.mesh();

    // Collect internal and boundary values
    List<List<Type>> stencilFld;
    collectData(fld, stencilFld);

    tmp<WeightedFieldType> twf
    (
        new WeightedFieldType
        (
            IOobject
            (
                fld.name(),
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensioned<WeightedType>
            (
                fld.name(),
                fld.dimensions(),
                Zero
            )
        )
    );
    WeightedFieldType& wf = twf();

    forAll(wf, celli)
    {
        const List<Type>& stField = stencilFld[celli];
        const List<WeightType>& stWeight = stencilWeights[celli];

        forAll(stField, i)
        {
            wf[celli] += stWeight[i]*stField[i];
        }
    }

    // Boundaries values?

    return twf;
}


// ************************************************************************* //
