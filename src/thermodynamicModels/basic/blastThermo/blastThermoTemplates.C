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

#include "blastThermo.H"
#include "wordIOList.H"
#include "compileTemplate.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, class Table>
typename Table::iterator Foam::blastThermo::lookupCstrIter
(
    const dictionary& thermoDict,
    Table* tablePtr,
    const int nCmpt,
    const char* cmptNames[],
    const word& thermoTypeName
)
{
    // Lookup the thermo package
    typename Table::iterator cstrIter = tablePtr->find(thermoTypeName);

    // Print error message if package not found in the table
    if (cstrIter == tablePtr->end())
    {
        const fileName origCODE_TEMPLATE_DIR(getEnv("FOAM_CODE_TEMPLATES"));
        fileName tempDir(Thermo::templateDir());
        setEnv("FOAM_CODE_TEMPLATES", tempDir, true);

        const word type(thermoDict.lookup("type"));
        word compileType = type + Thermo::typeName.capitalise() + "BlastThermo";

        Info<< "Trying to compile based on "
            << tempDir/compileType
            << endl;

        if
        (
            dynamicCode::allowSystemOperations
         && !dynamicCode::resolveTemplate(compileType).empty()
        )
        {
            dictionary thermoTypeDict;

            List<Pair<word>> entries;
            if (type == "detonating")
            {
                const dictionary& uDict =
                    thermoDict.subDict("reactants").subDict("thermoType");
                const dictionary& rDict =
                    thermoDict.subDict("products").subDict("thermoType");
                thermoTypeDict.add("reactants", uDict);
                thermoTypeDict.add("reactants", rDict);

                entries =
                    {
                        {"uTransport", uDict.lookup("transport")},
                        {"rTransport", rDict.lookup("transport")},
                        {"uThermo", uDict.lookup("thermo")},
                        {"rThermo", rDict.lookup("thermo")},
                        {"uEquationOfState", uDict.lookup("equationOfState")},
                        {"rEquationOfState", rDict.lookup("equationOfState")},
                        {"specie", "specieBlast"}
                    };
            }
            else
            {
                thermoTypeDict = thermoDict.subDict("thermoType");

                entries =
                    {
                        {"transport", thermoTypeDict.lookup("transport")},
                        {"thermo", thermoTypeDict.lookup("thermo")},
                        {"equationOfState", thermoTypeDict.lookup("equationOfState")},
                        {"specie", "specieBlast"}
                    };
            }
            compileTemplate thermo
            (
                compileType,
                thermoTypeName,
                entries
            );

            cstrIter = tablePtr->find(thermoTypeName);

            if (cstrIter == tablePtr->end())
            {
                FatalErrorInFunction
                    << "Compilation and linkage of "
                    << compileType << " type " << nl
                    << "thermoType" << thermoTypeDict << nl << nl
                    << "failed." << nl << nl
                    << "Valid " << Thermo::typeName << " types are:"
                    << nl << nl;
            }
        }
        else
        {
            // Print error message if package not found in the table
            FatalErrorInFunction
                << "Unknown " << Thermo::typeName << " type " << nl
                << "thermoType " << thermoTypeName << nl << nl
                << "Valid " << Thermo::typeName << " types are:"
                << nl << nl;
        }

        if (!origCODE_TEMPLATE_DIR.empty())
        {
            system
            (
                word
                (
                    "export FOAM_CODE_TEMPLATES=" + origCODE_TEMPLATE_DIR
                ).c_str()
            );
        }

        if (cstrIter == tablePtr->end())
        {
            // Get the list of all the suitable thermo packages available
            wordList validThermoTypeNames
            (
                tablePtr->sortedToc()
            );

            // Build a table of the thermo packages constituent parts
            // Note: row-0 contains the names of constituent parts
            List<wordList> validThermoTypeNameCmpts
            (
                validThermoTypeNames.size() + 1
            );

            validThermoTypeNameCmpts[0].setSize(6);
            forAll(validThermoTypeNameCmpts[0], j)
            {
                validThermoTypeNameCmpts[0][j] = cmptNames[j];
            }

            // Split the thermo package names into their constituent parts
            // Removing incompatible entries from the list
            label j = 1;
            forAll(validThermoTypeNames, i)
            {
                wordList names
                (
                    blastThermo::splitThermoName(validThermoTypeNames[i])
                );

                if (names.size())
                {
                    validThermoTypeNameCmpts[j++] = names;
                }
            }
            validThermoTypeNameCmpts.setSize(j);

            // Print the table of available packages
            // in terms of their constituent parts
            printTable(validThermoTypeNameCmpts, FatalError);

            FatalError<< exit(FatalError);
        }
    }

    return cstrIter;
}


template<class Thermo, class Table>
typename Table::iterator Foam::blastThermo::lookupCstrIter
(
    const dictionary& thermoDict,
    Table* tablePtr
)
{
    word state
    (
        thermoDict.lookupOrDefault<word>("state", Thermo::typeName)
    );
    word type(thermoDict.lookup<word>("type"));
    word typeName(state + "<" + type + "<");
    word thermoTypeName;
    if (type == "detonating")
    {
        const dictionary& uThermoTypeDict
        (
            thermoDict.subDict("reactants").subDict("thermoType")
        );
        const dictionary& rThermoTypeDict
        (
            thermoDict.subDict("products").subDict("thermoType")
        );

        const word uthermoTypeName(readThermoType(uThermoTypeDict));
        const word rthermoTypeName(readThermoType(rThermoTypeDict));
        thermoTypeName = uthermoTypeName + ',' + rthermoTypeName;
    }
    else
    {

        const dictionary& thermoTypeDict(thermoDict.subDict("thermoType"));
        thermoTypeName = readThermoType(thermoTypeDict);
    }

    typeName += thermoTypeName + ">>";
    Info<< "Selecting thermodynamics package " << typeName << endl;

    const int nCmpt = 6;
    const char* cmptNames[nCmpt] =
    {
        "state",
        "type",
        "transport",
        "thermo",
        "equationOfState",
        "specie"
    };

    return lookupCstrIter<Thermo, Table>
    (
        thermoDict,
        tablePtr,
        nCmpt,
        cmptNames,
        typeName
    );
}


template<class Thermo>
Foam::autoPtr<Thermo> Foam::blastThermo::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
{
    typename Thermo::dictionaryConstructorTable::iterator cstrIter =
        lookupCstrIter<Thermo, typename Thermo::dictionaryConstructorTable>
        (
            dict,
            Thermo::dictionaryConstructorTablePtr_
        );

    return autoPtr<Thermo>(cstrIter()(mesh, dict, phaseName, masterName));
}

// ************************************************************************* //
