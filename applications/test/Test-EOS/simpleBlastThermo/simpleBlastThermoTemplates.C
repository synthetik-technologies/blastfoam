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

#include "simpleBlastThermo.H"
#include "compileTemplate.H"
#include "OSspecific.H"

/* * * * * * * * * * * * * * * * * * static data * * * * * * * * * * * * * * */

template<class Thermo, class Table>
typename Table::iterator Foam::simpleBlastThermo::lookupCstrIter
(
    const dictionary& thermoTypeDict,
    Table* tablePtr,
    const int nCmpt,
    const char* cmptNames[],
    const word& thermoTypeName
)
{
    // Lookup the thermo package
    typename Table::iterator cstrIter =
        tablePtr->find(thermoTypeName);

        // Print error message if package not found in the table
    if (cstrIter == tablePtr->end())
    {
        const word compileType("simpleBlastThermo");
        const fileName origCODE_TEMPLATE_DIR(getEnv("FOAM_CODE_TEMPLATES"));
        const fileName BLAST_DIR(getEnv("BLAST_DIR"));
        setEnv("FOAM_CODE_TEMPLATES", BLAST_DIR/"etc/codeTemplates", true);

        if
        (
            dynamicCode::allowSystemOperations
         && !dynamicCode::resolveTemplate(compileType).empty()
        )
        {
            compileTemplate thermo
            (
                compileType,
                thermoTypeName,
                List<Pair<word>>
                (
                    {
                        {"transport", thermoTypeDict.lookup("transport")},
                        {"thermo", thermoTypeDict.lookup("thermo")},
                        {
                            "equationOfState",
                            thermoTypeDict.lookup("equationOfState")
                        },
                        {"specie", "specieBlast"}
                    }
                )
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
            FatalError
                << "Unknown " << Thermo::typeName << " type " << nl
                << "thermoType " << thermoTypeName << nl << nl
                << "Valid " << Thermo::typeName << " types are:"
                << nl << nl;

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

            validThermoTypeNameCmpts[0].setSize(nCmpt);
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
                    splitThermoName(validThermoTypeNames[i])
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
typename Table::iterator Foam::simpleBlastThermo::lookupCstrIter
(
    const dictionary& thermoDict,
    Table* tablePtr
)
{
    const dictionary& thermoTypeDict(thermoDict.subDict("thermoType"));
    word thermoTypeName = readThermoType(thermoTypeDict);

    Info<< "Selecting thermodynamics package " << thermoTypeName << endl;

    const int nCmpt = 4;
    const char* cmptNames[nCmpt] =
    {
        "transport",
        "thermo",
        "equationOfState",
        "specie"
    };

    return lookupCstrIter<Thermo, Table>
    (
        thermoTypeDict,
        tablePtr,
        nCmpt,
        cmptNames,
        thermoTypeName
    );
}

// ************************************************************************* //
