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
    Template for use with dynamic code generation of a Function3.

    - without state

SourceFiles
    Function3Template.C

\*---------------------------------------------------------------------------*/

#ifndef Function3Template_H
#define Function3Template_H

#include "Function3.H"

//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function3s
{

/*---------------------------------------------------------------------------*\
                           A templated Function3
\*---------------------------------------------------------------------------*/

class ${typeName}Function3${TemplateType}
:
    public FieldFunction3<${TemplateType}, ${typeName}Function3${TemplateType}>
{

public:

    // Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from entry name and dictionary
        ${typeName}Function3${TemplateType}
        (
            const word& entryName,
            const dictionary& dict
        );

        //- Copy constructor
        ${typeName}Function3${TemplateType}
        (
            const ${typeName}Function3${TemplateType}& f1
        );

        //- Construct and return a clone
        virtual tmp<Function3<${TemplateType}>> clone() const
        {
            return tmp<Function3<${TemplateType}>>
            (
                new ${typeName}Function3${TemplateType}(*this)
            );
        }


    //- Destructor
    virtual ~${typeName}Function3${TemplateType}();


    // Member Functions

        //- Return value as a function of two scalars
        inline virtual ${TemplateType} value
        (
            const scalar x,
            const scalar y
        ) const
        {
//{{{ begin code
    ${code}
//}}} end code
        }

        //- Write data to dictionary stream
        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ${typeName}Function3${TemplateType}&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function3s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
