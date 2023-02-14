/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

#include "PatchFieldSelector.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::PtrList<Foam::PatchFieldSelector<Type>>
Foam::PatchFieldSelector<Type>::New
(
    const fvPatchField<Type>& pfield,
    const dictionary& dict
)
{
    PtrList<PatchFieldSelector<Type>> pacceleration;
    if
    (
        isA<calculatedFvPatchField<Type>>(pfield)
     || isA<fixedValueFvPatchField<Type>>(pfield)
    )
    {
        DebugInfo
            << "Using value type selector for patch "
            << string(pfield.patch().name())
            << " and field " << string(pfield.internalField().name())
            << endl;
        pacceleration.setSize(1);
        pacceleration.set
        (
            0,
            new ValueFvPatchFieldSelector<Type>(pfield)
        );
    }
    else if (isA<fixedGradientFvPatchField<Type>>(pfield))
    {
        DebugInfo
            << "Using gradient type selector for patch "
            << string(pfield.patch().name())
            << " and field " << string(pfield.internalField().name())
            << endl;
        pacceleration.setSize(1);
        pacceleration.set
        (
            0,
            new FixedGradientFvPatchFieldSelector<Type>(pfield)
        );
    }
    else if (isA<mixedFvPatchField<Type>>(pfield))
    {
        bool useRefValue = true;
        bool useRefGrad = true;
        bool useValueFraction = false;
        if (dict.isDict(pfield.patch().name()))
        {
            const dictionary& patchDict = dict.subDict(pfield.patch().name());
            useRefValue =
                patchDict.template lookupOrDefault("useRefValue", true);
            useRefGrad =
                patchDict.template lookupOrDefault("useRefGrad", true);
            // useValueFraction =
            //     patchDict.template lookupOrDefault("useValueFraction", true);
        }
        label nFields =
            label(useRefValue) + label(useRefGrad) + label(useValueFraction);

        pacceleration.setSize(nFields);
        label selectori = 0;
        if (useRefValue)
        {
            DebugInfo
                << "Using mixed type selector with a refValue for patch "
                << string(pfield.patch().name())
                << " and field " << string(pfield.internalField().name())
                << endl;
            pacceleration.set
            (
                selectori++,
                new MixedFvPatchFieldSelector<Type>
                (
                    pfield,
                    MixedFvPatchFieldSelector<Type>::REF_VALUE
                )
            );
        }
        if (useRefGrad)
        {
            DebugInfo
                << "Using mixed type selector with a refGrad for patch "
                << string(pfield.patch().name())
                << " and field " << string(pfield.internalField().name())
                << endl;
            pacceleration.set
            (
                selectori++,
                new MixedFvPatchFieldSelector<Type>
                (
                    pfield,
                    MixedFvPatchFieldSelector<Type>::REF_GRADIENT
                )
            );
        }
        // if (useValueFraction)
        // {
        //     pacceleration.set
        //     (
        //         selectori++,
        //         new MixedFvPatchFieldSelector<Type>
        //         (
        //             pfield,
        //             MixedFvPatchFieldSelector<Type>::VALUE_FRACTION
        //         )
        //     );
        // }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported fvPatchField type for "
            << pfield.internalField().name() << " of type "
            << pfield.type()
            << ". The following base types are supported: " << nl
            << "    calculated" << nl
            << "    fixedValue" << nl
            << "    fixedGradient" << nl
            << "    mixed" << nl
            << abort(FatalError);
    }
    return pacceleration;
}


template<class Type>
Foam::PtrList<Foam::PatchFieldSelector<Type>>
Foam::PatchFieldSelector<Type>::New
(
    const fvsPatchField<Type>& pfield,
    const dictionary& dict
)
{
    PtrList<PatchFieldSelector<Type>> pacceleration;
    if
    (
        isA<calculatedFvsPatchField<Type>>(pfield)
     || isA<fixedValueFvsPatchField<Type>>(pfield)
    )
    {
        DebugInfo
            << "Using value type selector for patch "
            << string(pfield.patch().name())
            << " and field " << string(pfield.internalField().name())
            << endl;
        pacceleration.setSize(1);
        pacceleration.set
        (
            0,
            new ValueFvsPatchFieldSelector<Type>(pfield)
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported fvsPatchField type for "
            << pfield.internalField().name() << " of type "
            << pfield.type()
            << ". The following base types are supported: " << nl
            << "    calculated" << nl
            << "    fixedValue" << nl
            << abort(FatalError);
    }
    return pacceleration;
}



template<class Type>
Foam::PtrList<Foam::PatchFieldSelector<Type>>
Foam::PatchFieldSelector<Type>::New
(
    const pointPatchField<Type>& pfield,
    const dictionary& dict
)
{
    PtrList<PatchFieldSelector<Type>> pacceleration;
    if (isA<valuePointPatchField<Type>>(pfield))
    {
        if (debug)
        {
            const polyPatch& pp =
                pfield.internalField().mesh().mesh().boundaryMesh()
                [
                    pfield.patch().index()
                ];
            Info<< "Using value type selector for patch "
                << string(pp.name())
                << " and field " << string(pfield.internalField().name())
                << endl;
        }
        pacceleration.setSize(1);
        pacceleration.set
        (
            0,
            new ValuePointPatchFieldSelector<Type>(pfield)
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported pointPatchField type for "
            << pfield.internalField().name() << " of type "
            << pfield.type()
            << ". The following base types are supported: " << nl
            << "    value" << nl
            << "    fixedValue" << nl
            << abort(FatalError);
    }
    return pacceleration;
}


// ************************************************************************* //
