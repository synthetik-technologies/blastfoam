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

#include "basicThermoModel.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Table>
typename Table::iterator Foam::basicThermoModel::lookupThermo
(
    Table* tablePtr,
    const char* cmptNames[],
    const word& thermoTypeName
)
{
    // Lookup the thermo package
    typename Table::iterator cstrIter = tablePtr->find(thermoTypeName);

    // Print error message if package not found in the table
    if (cstrIter == tablePtr->end())
    {
        FatalErrorInFunction
            << "Unknown " << basicThermoModel::typeName << " type " << nl
            << "thermoType " << thermoTypeName << nl << nl
            << "Valid " << basicThermoModel::typeName << " types are:"
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

        validThermoTypeNameCmpts[0].setSize(5);
        forAll(validThermoTypeNameCmpts[0], j)
        {
            validThermoTypeNameCmpts[0][j] = cmptNames[j];
        }

        // Split the thermo package names into their constituent parts
        // Removing incompatible entries from the list
        label j = 0;
        forAll(validThermoTypeNames, i)
        {
            wordList names
            (
                basicThermoModel::splitThermoName(validThermoTypeNames[i])
            );

            if (names.size())
            {
                validThermoTypeNameCmpts[++j] = names;
            }
        }
        validThermoTypeNameCmpts.setSize(j);

        // Print the table of available packages
        // in terms of their constituent parts
        printTable(validThermoTypeNameCmpts, FatalError);

        FatalError<< exit(FatalError);
    }

    return cstrIter;
}


template<class Table>
typename Table::iterator Foam::basicThermoModel::lookupThermo
(
    const dictionary& thermoDict,
    Table* tablePtr
)
{
    word thermoTypeName;
    Switch detonating(thermoDict.lookupType<word>("type") == "detonating");
    if (detonating)
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
        thermoTypeName =
            word("detonating<")
          + uthermoTypeName
          + ','
          + rthermoTypeName
          + '>';
    }
    else
    {

        const dictionary& thermoTypeDict(thermoDict.subDict("thermoType"));
        thermoTypeName = readThermoType(thermoTypeDict);
    }
    Info<< "Selecting thermodynamics package " << thermoTypeName << endl;

    const int nCmpt = 5;
    const char* cmptNames[nCmpt] =
    {
        "detonating",
        "transport",
        "thermo",
        "equationOfState",
        "specie"
    };

    return lookupThermo<Table>
    (
        tablePtr,
        cmptNames,
        thermoTypeName
    );
}
// ************************************************************************* //
