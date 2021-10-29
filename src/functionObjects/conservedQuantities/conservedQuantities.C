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

    const volScalarField::Internal& V(this->mesh_.V());
    forAll(fields_, i)
    {
        const word& name = fields_[i];
        if (obr_.foundObject<volScalarField>(name))
        {
            const volScalarField& f = obr_.lookupObject<volScalarField>(name);
            scalars_.insert(name, (sum(f()*V)));
            scalar0s_.insert(name, (sum(f()*V)));
        }
        else if (obr_.foundObject<volVectorField>(name))
        {
            const volVectorField& f = obr_.lookupObject<volVectorField>(name);
            vectors_.insert(name, (sum(f()*V)));
            vector0s_.insert(name, (sum(f()*V)));
        }
        else if (obr_.foundObject<volSymmTensorField>(name) )
        {
            const volSymmTensorField& f =
                obr_.lookupObject<volSymmTensorField>(name);
            symmTensors_.insert(name, (sum(f()*V)));
            symmTensor0s_.insert(name, (sum(f()*V)));
        }
        else if (obr_.foundObject<volSphericalTensorField>(name) )
        {
            const volSphericalTensorField& f =
                obr_.lookupObject<volSphericalTensorField>(name);
            sphTensors_.insert
            (
                name,
                (sum(f()*V))
            );
            sphTensor0s_.insert
            (
                name,
                sum(f()*V)
            );
        }
        else if (obr_.foundObject<volTensorField>(name) )
        {
            const volTensorField& f =
                obr_.lookupObject<volTensorField>(name);
            tensors_.insert(name, (sum(f()*V)));
            tensor0s_.insert(name, (sum(f()*V)));
        }
        else
        {
            WarningInFunction
                << name << " was not found. Skipping" << endl;
        }
    }

    Log << endl;

    return true;
}


bool Foam::functionObjects::conservedQuantities::execute()
{
    const volScalarField::Internal& V(this->mesh_.V());
    forAll(fields_, i)
    {
        const word& name = fields_[i];
        if (obr_.foundObject<volScalarField>(name))
        {
            const volScalarField& f = obr_.lookupObject<volScalarField>(name);
            scalars_[name] = sum(f()*V);
            scalar s = scalars_[name].value();
            scalar s0 = scalar0s_[name].value();
            Info << name << " " << f.dimensions() << " :" << nl
                <<"    Original = " << s0 << nl
                <<"    Current = " << s << nl
                <<"    difference (abs/rel) = " << (s - s0) << ", "
                << mag(s - s0)/s0
                << endl;
        }
        else if (obr_.foundObject<volVectorField>(name))
        {
            const volVectorField& f = obr_.lookupObject<volVectorField>(name);
            vectors_[name] = sum(f()*V);
            vector v = vectors_[name].value();
            vector v0 = vector0s_[name].value();
            Info << name << " " << f.dimensions() << " :" << nl
                <<"    Original = " << v0 << nl
                <<"    Current = " << v << nl
                <<"    difference (abs/rel) = "<< (v - v0) << ", "
                << mag(v - v0)/mag(v0)
                << endl;
        }
        else if (obr_.foundObject<volSymmTensorField>(name) )
        {
            const volSymmTensorField& f =
                obr_.lookupObject<volSymmTensorField>(name);
            symmTensors_[name] = sum(f()*V);
            symmTensor st = symmTensors_[name].value();
            symmTensor st0 = symmTensor0s_[name].value();
            Info << name << " " << f.dimensions() << " :" << nl
                <<"    Original = " << st0 << nl
                <<"    Current = " << st << nl
                <<"    difference (abs/rel) = " << (st - st0) << ", "
                << mag(st - st0)/mag(st0)
                << endl;
        }
        else if (obr_.foundObject<volSphericalTensorField>(name) )
        {
            const volSphericalTensorField& f =
                obr_.lookupObject<volSphericalTensorField>(name);
            sphTensors_[name] = sum(f()*V);
            symmTensor st = sphTensors_[name].value();
            symmTensor st0 = sphTensor0s_[name].value();
            Info << name << " " << f.dimensions() << " :" << nl
                <<"    Original = " << st0 << nl
                <<"    Current = " << st << nl
                <<"    difference (abs/rel) = " << (st - st0) << ", "
                << mag(st - st0)/mag(st0)
                << endl;
        }
        else if (obr_.foundObject<volTensorField>(name) )
        {
            const volTensorField& f =
                obr_.lookupObject<volTensorField>(name);
            tensors_[name] = sum(f()*V);
            tensor t = tensors_[name].value();
            tensor t0 = tensor0s_[name].value();
            Info << name << " " << f.dimensions() << " :" << nl
                <<"    Original = " << t0 << nl
                <<"    Current = " << t << nl
                <<"    difference (abs/rel) = " << (t - t0) << ", "
                << mag(t - t0)/mag(t0)
                << endl;
        }
    }
    return true;
}


bool Foam::functionObjects::conservedQuantities::write()
{
    return true;
}


// ************************************************************************* //
