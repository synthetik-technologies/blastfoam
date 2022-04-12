/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "CodedMultivariateEquationTemplate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        ${typeName}_${TemplateType}MultivariateEquation,
        0
    );
    regEquation<${TemplateType}, MultivariateEquation>::
    adddictionaryConstructorToTable<${typeName}_${TemplateType}MultivariateEquation>
        ${typeName}_${TemplateType}RegMultivariateEquationConstructorToTable_;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ${SHA1sum}
    //
    // Unique function name that can be checked if the correct library version
    // has been loaded
    void ${typeName}_${SHA1sum}(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::${typeName}_${TemplateType}MultivariateEquation::
${typeName}_${TemplateType}MultivariateEquation
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    regEquation<${TemplateType}, MultivariateEquation>(obr, dict)
{
    if (${verbose:-false})
    {
        Info<< "Construct ${typeName} sha1: ${SHA1sum} from dictionary\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::${typeName}_${TemplateType}MultivariateEquation::
~${typeName}_${TemplateType}MultivariateEquation()
{
    if (${verbose:-false})
    {
        Info<< "Destroy ${typeName} sha1: ${SHA1sum}\n";
    }
}

// ************************************************************************* i/
