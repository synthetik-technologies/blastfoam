/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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

#include "basicBlastThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(basicBlastThermo, 0);
    defineRunTimeSelectionTable(basicBlastThermo, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicBlastThermo::basicBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    timeIntegrationSystem
    (
        IOobject::groupName("basicBlastThermo", phaseName),
        mesh
    ),
    masterName_(masterName),
    phaseName_(phaseName),
    rho_
    (
        mesh.lookupObjectRef<volScalarField>
        (
            IOobject::groupName("rho", phaseName)
        )
    ),
    T_(mesh.lookupObject<volScalarField>(IOobject::groupName("T", masterName_))),
    e_(mesh.lookupObject<volScalarField>(IOobject::groupName("e", masterName_))),
    residualAlpha_("residualAlpha", dimless, 0.0),
    residualRho_("residualRho", dimDensity, 0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicBlastThermo::~basicBlastThermo()
{}

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicBlastThermo> Foam::basicBlastThermo::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
{
    dictionaryConstructorTable::iterator cstrIter =
        lookupCstrIter<basicBlastThermo, dictionaryConstructorTable>
        (
            dict,
            dictionaryConstructorTablePtr_
        );

    return cstrIter()
    (
        mesh,
        dict,
        phaseName,
        masterName
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::UIndirectList<Foam::scalar> Foam::basicBlastThermo::cellSetScalarList
(
    const volScalarField& psi,
    const labelList& cells
)
{
    return UIndirectList<scalar>(psi, cells);
}


Foam::word Foam::basicBlastThermo::readThermoType(const dictionary& dict)
{
    return word
    (
        word(dict.lookup("transport")) + '<'
      + word(dict.lookup("thermo")) + '<'
      + word(dict.lookup("equationOfState"))  + '<'
      + word("specie") + ">>>"
    );
}

Foam::wordList Foam::basicBlastThermo::splitThermoName
(
    const word& thermoName,
    const int nCmpt
)
{
    wordList cmpts(nCmpt);

    string::size_type beg=0, end=0, endb=0, endc=0;
    int i = 0;

    while
    (
        (endb = thermoName.find('<', beg)) != string::npos
     || (endc = thermoName.find(',', beg)) != string::npos
    )
    {
        if (endb == string::npos)
        {
            end = endc;
        }
        else if ((endc = thermoName.find(',', beg)) != string::npos)
        {
            end = std::min(endb, endc);
        }
        else
        {
            end = endb;
        }

        if (beg < end)
        {
            cmpts[i] = thermoName.substr(beg, end-beg);
            cmpts[i++].replaceAll(">","");

            // If the number of number of components in the name
            // is greater than nCmpt return an empty list
            if (i == nCmpt)
            {
                return wordList();
            }
        }
        beg = end + 1;
    }

    // If the number of number of components in the name is not equal to nCmpt
    // return an empty list
    if (i + 1 != nCmpt)
    {
        return wordList();
    }

    if (beg < thermoName.size())
    {
        cmpts[i] = thermoName.substr(beg, string::npos);
        cmpts[i].replaceAll(">","");
    }

    return cmpts;
}


Foam::wordList Foam::basicBlastThermo::splitThermoName
(
    const word& thermoName
)
{
    wordList cmpts;
    string::size_type beg=0, end=0, endb=0, endc=0;

    while
    (
        (endb = thermoName.find('<', beg)) != string::npos
     || (endc = thermoName.find(',', beg)) != string::npos
    )
    {
        if (endb == string::npos)
        {
            end = endc;
        }
        else if ((endc = thermoName.find(',', beg)) != string::npos)
        {
          end = std::min(endb, endc);
        }
        else
        {
            end = endb;
        }

        if (beg < end)
        {
            word newStr = thermoName.substr(beg, end-beg);
            newStr.replaceAll(">","");
            cmpts.append(newStr);
        }
        beg = end + 1;
    }
    if (beg < thermoName.size())
    {
        word newStr = thermoName.substr(beg, string::npos);
        newStr.replaceAll(">","");
        cmpts.append(newStr);
    }

    wordList cmptsFinal(6);
    if (cmpts[0] == "detonating")
    {
        cmptsFinal[0] = cmpts[0];
        cmptsFinal[1] = cmpts[1];
        cmptsFinal[2] = cmpts[2] + '/' + cmpts[6];
        cmptsFinal[3] = cmpts[3] + '/' + cmpts[7];
        cmptsFinal[4] = cmpts[4] + '/' + cmpts[8];
        cmptsFinal[5] = cmpts[5] + '/' + cmpts[9];
    }
    else
    {
        cmptsFinal[0] = cmpts[0];
        cmptsFinal[1] = cmpts[1];
        cmptsFinal[2] = cmpts[2];
        cmptsFinal[3] = cmpts[3];
        cmptsFinal[4] = cmpts[4];
        cmptsFinal[5] = cmpts[5];
    }

    return cmptsFinal;
}


Foam::volScalarField& Foam::basicBlastThermo::lookupOrConstruct
(
    const fvMesh& mesh,
    const word& name,
    const IOobject::readOption rOpt,
    const IOobject::writeOption wOpt,
    const dimensionSet& dims,
    const bool allowNoGroup
)
{
    const word baseName(IOobject::member(name));
    if (!mesh.foundObject<volScalarField>(name))
    {
        volScalarField* fPtr = nullptr;
        IOobject io
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            wOpt
        );
        IOobject baseIo
        (
            baseName,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            wOpt
        );

        if (io.typeHeaderOk<volScalarField>(true))
        {
            io.readOpt() = IOobject::MUST_READ;
            fPtr =
                new volScalarField
                (
                    io,
                    mesh
                );
        }
        else if (mesh.foundObject<volScalarField>(baseName) && allowNoGroup)
        {
            const volScalarField& baseField =
                mesh.lookupObjectRef<volScalarField>(baseName);
            fPtr =
                new volScalarField
                (
                    io,
                    baseField,
                    baseField.boundaryField()
                );
        }
        else if (baseIo.typeHeaderOk<volScalarField>(true) && allowNoGroup)
        {
            baseIo.readOpt() = IOobject::MUST_READ;
            fPtr =
                new volScalarField
                (
                    baseIo,
                    mesh
                );

            // Rename to the desired name
            fPtr->rename(name);
        }
        else if (rOpt != IOobject::MUST_READ)
        {
            fPtr =
                new volScalarField
                (
                    io,
                    mesh,
                    dimensionedScalar("0", dims, Zero)
                );
        }
        else if (allowNoGroup)
        {
            FatalErrorInFunction
                << name << " or " << baseName << " is required" << endl
                << abort(FatalError);
        }
        else
        {
            FatalErrorInFunction
                << name << " is required" << endl
                << abort(FatalError);
        }

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return mesh.lookupObjectRef<volScalarField>(name);
}


// ************************************************************************* //
