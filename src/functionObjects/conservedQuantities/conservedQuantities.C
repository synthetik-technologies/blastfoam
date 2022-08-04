/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
28-10-2020 Synthetik Applied Technologies: | Calculate conservedQuantities
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "conservedQuantities.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(conservedQuantities, 0);
    addToRunTimeSelectionTable(functionObject, conservedQuantities, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::conservedQuantities::conservedQuantities
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fields_(dict.lookup("fields"))
{
    read(dict);

    const volScalarField::Internal& V(this->mesh_.V());
    forAll(fields_, i)
    {
        const word& name = fields_[i];

        bool found = false;
        #define CreateConservedVariables(Type, none)                    \
        if                                                              \
        (                                                               \
            obr_.foundObject                                            \
            <                                                           \
                GeometricField<Type, fvPatchField, volMesh>             \
            >(name)                                                     \
        )                                                               \
        {                                                               \
            const GeometricField<Type, fvPatchField, volMesh>& f =      \
                obr_.lookupObject                                       \
                <                                                       \
                    GeometricField<Type, fvPatchField, volMesh>         \
                >(name);                                                \
            dimensioned<Type> fV(sum(f()*V));                           \
            vals0<Type>().insert(name, fV);                             \
            vals<Type>().insert(name, fV);                              \
            found = true;                                               \
        }

        FOR_ALL_FIELD_TYPES(CreateConservedVariables);

        #undef CreateConservedVariables

        if (!found)
        {
            cannotFindObject(name);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::conservedQuantities::~conservedQuantities()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::conservedQuantities::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Log << type() << " " << this->name() << ":" << nl;

    fields_ = wordList(dict.lookup("fields"));

    Log << endl;

    return true;
}


bool Foam::functionObjects::conservedQuantities::execute()
{
    const volScalarField::Internal& V(this->mesh_.V());
    forAll(fields_, i)
    {
        const word& name = fields_[i];

        bool found = false;
        #define UpdateConservedVariables(Type, none)                    \
        if                                                              \
        (                                                               \
            obr_.foundObject                                            \
            <                                                           \
                GeometricField<Type, fvPatchField, volMesh>             \
            >(name)                                                     \
        )                                                               \
        {                                                               \
            const GeometricField<Type, fvPatchField, volMesh>& f =      \
                obr_.lookupObject                                       \
                <                                                       \
                    GeometricField<Type, fvPatchField, volMesh>         \
                >(name);                                                \
                                                                        \
            dimensioned<Type>& fV = vals<Type>()[name];                 \
            fV = sum(f()*V);                                            \
            const Type& val0 = vals0<Type>()[name].value();             \
            const Type& val = fV.value();                               \
                                                                        \
            Info<< name << " " << f.dimensions() << " :" << nl          \
                << "    Original = " << val0 << nl                      \
                << "    Current = " << val << nl                        \
                << "    difference (abs/rel) = " << (val - val0)        \
                << ", "  << mag(val - val0)/max(mag(val0), small)       \
                << endl;                                                \
                                                                        \
            found = true;                                               \
        }

        FOR_ALL_FIELD_TYPES(UpdateConservedVariables);

        #undef UpdateConservedVariables

        if (!found)
        {
            cannotFindObject(name);
        }
    }
    return true;
}


bool Foam::functionObjects::conservedQuantities::write()
{
    return true;
}


// ************************************************************************* //
