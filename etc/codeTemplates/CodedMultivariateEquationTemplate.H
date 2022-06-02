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

Description
    Template for use with dynamic code generation of an Equation

SourceFiles
    CodedMultivariateTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef CodedMultivariateTemplate_H
#define CodedMultivariateTemplate_H

#include "MultivariateEquationsFwd.H"
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

class ${typeName}_${TemplateType}MultivariateEquation
:
    public regEquation<${TemplateType}, MultivariateEquation>
{
public:

    // Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from entry name and dictionary
        ${typeName}_${TemplateType}MultivariateEquation
        (
            const objectRegistry&,
            const dictionary&
        );



    //- Destructor
    virtual ~${typeName}_${TemplateType}MultivariateEquation();


    // Member Functions

        //- Return the number of derivatives
        inline virtual label nDerivatives() const
        {
            return ${nDerivatives};
        }

        //- Return value as a function of scalar variable
        inline virtual void FX
        (
            const UList<scalar>& x,
            const label li,
            List<${TemplateType}>& fx
        ) const
        {
//{{{ begin code
    ${fx_code}
//}}} end code
        }

        //- Return value as a function of scalar variable
        inline virtual void jacobian
        (
            const UList<scalar>& x,
            const label li,
            List<${TemplateType}>& fx,
            RectangularMatrix<${TemplateType}>& J
        ) const
        {
//{{{ begin code
    ${dfdx_code}
//}}} end code
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
