/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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

Description
    Template for use with dynamic code generation of an Equation

SourceFiles
    CodedEquationTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef CodedEquationTemplate_H
#define CodedEquationTemplate_H

#include "EquationsFwd.H"
#include "regEquation.H"

//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           A templated Equation
\*---------------------------------------------------------------------------*/

class ${typeName}_${TemplateType}Equation
:
    public regEquation<${TemplateType}, Equation>
{
public:

    // Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from entry name and dictionary
        ${typeName}_${TemplateType}Equation
        (
            const objectRegistry& obr,
            const dictionary&
        );



    //- Destructor
    virtual ~${typeName}_${TemplateType}Equation();


    // Member Functions

        //- Return the number of derivatives
        inline virtual label nDerivatives() const
        {
            return ${nDerivatives};
        }
        //- Return value as a function of scalar variable
        inline virtual ${TemplateType} fx(const scalar x, const label li) const
        {
//{{{ begin code
    ${fx_code}
//}}} end code
        }

        //- Return value as a function of scalar variable
        inline virtual ${TemplateType} dfdx(const scalar x, const label li) const
        {
//{{{ begin code
    ${dfdx_code}
//}}} end code
        }

        //- Return value as a function of scalar variable
        inline virtual ${TemplateType} d2fdx2(const scalar x, const label li) const
        {
//{{{ begin code
    ${d2fdx2_code}
//}}} end code
        }

        //- Return value as a function of scalar variable
        inline virtual ${TemplateType} d3fdx3(const scalar x, const label li) const
        {
//{{{ begin code
    ${d3fdx3_code}
//}}} end code
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
