/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const solidModel& lookupSolidModel(const objectRegistry& obr)
{
    return lookupSolidModel(obr, obr.name());
}


const DSolidModel& lookupDSolidModel(const objectRegistry& obr)
{
    return lookupDSolidModel(obr, obr.name());
}

const solidModel& lookupSolidModel
(
    const objectRegistry& obr,
    const word& baseMeshRegionName
)
{
    if (obr.foundObject<solidModel>("solidProperties"))
    {
        return obr.lookupObject<solidModel>
        (
            "solidProperties"
        );
    }
    else if (obr.parent().foundObject<solidModel>("solidProperties"))
    {
        return obr.parent().lookupObject<solidModel>
        (
            "solidProperties"
        );
    }
    else
    {
        HashTable<const objectRegistry*> obrs
        (
            obr.parent().lookupClass<objectRegistry>()
        );

        forAllConstIter
        (
            HashTable<const objectRegistry*>,
            obrs,
            iter
        )
        {
            if
            (
                iter.key() == baseMeshRegionName
             && iter()->foundObject<solidModel>("solidProperties")
            )
            {
                return iter()->lookupObject<solidModel>
                (
                    "solidProperties"
                );
            }
        }
    }

    FatalErrorInFunction
        << "Could not find " << solidModel::typeName
        << "for region " << obr.name() << nl << nl
        << "solidModels in the objectRegistry: "
        << obr.lookupClass<solidModel>().toc() << nl << nl
        << "solidModels in the parent objectRegistry:"
        << obr.parent().lookupClass<solidModel>().toc() << abort(FatalError);

    // Keep the compiler happy
    return obr.lookupObject<solidModel>("none");
}


const DSolidModel& lookupDSolidModel
(
    const objectRegistry& obr,
    const word& baseMeshRegionName
)
{
    return refCast<const DSolidModel>(lookupSolidModel(obr, baseMeshRegionName));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
