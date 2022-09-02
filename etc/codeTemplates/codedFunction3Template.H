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
    Template for use with dynamic code generation of a Function3.

    - without state

SourceFiles
    Function3Template.C

\*---------------------------------------------------------------------------*/

#ifndef codedFunction3Template_H
#define codedFunction3Template_H

#include "Function3.H"
#include "fieldTypes.H"

// Include all function1s for easier access
#include "None.H"
#include "Constant.H"
#include "Uniform.H"
#include "ZeroConstant.H"
#include "OneConstant.H"
#include "Polynomial1.H"
#include "Sine.H"
#include "Square.H"
#include "Table.H"
#include "UniformTable1.H"
#include "NonUniformTable1.H"
#include "Scale.H"
#include "CodedFunction1.H"

// Include all function2s for easier access
#include "None2.H"
#include "Constant2.H"
#include "ZeroConstant2.H"
#include "OneConstant2.H"
#include "Scale2.H"
#include "UniformTable2.H"
#include "CodedFunction2.H"

//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function3s
{

/*---------------------------------------------------------------------------*\
                    Helper functions to construct a Function1
\*---------------------------------------------------------------------------*/

// Helper functions to make using exisiting functions easier
#define Func10(Name)                                    \
    Function1<${TemplateType}>::New                     \
    (                                                   \
        #Name,                                          \
        dict_                                           \
    )()

#define Func11(Name, Func)                              \
    autoPtr<Function1<${TemplateType}>> Func##Ptr       \
    (                                                   \
        Function1<${TemplateType}>::New                 \
        (                                               \
            #Name,                                      \
            dict_                                       \
        )                                               \
    );                                                  \
    Function1<${TemplateType}>& Func = Func##Ptr();

#define GET_FUNC1_MACRO(_0,_1,NAME,...) NAME
#define Func1(...) \
    GET_FUNC1_MACRO(__VA_ARGS__, Func11, Func10)(__VA_ARGS__)


/*---------------------------------------------------------------------------*\
                    Helper functions to construct a Function2
\*---------------------------------------------------------------------------*/

#define Func20(Name)                                    \
    Function2<${TemplateType}>::New                     \
    (                                                   \
        #Name,                                          \
        dict_                                           \
    )()

#define Func21(Name, Func)                              \
    autoPtr<Function2<${TemplateType}>> Func##Ptr       \
    (                                                   \
        Function2<${TemplateType}>::New                 \
        (                                               \
            #Name,                                      \
            dict_                                       \
        )                                               \
    );                                                  \
    Function2<${TemplateType}>& Func = Func##Ptr();

#define GET_FUNC2_MACRO(_0,_1,NAME,...) NAME
#define Func2(...) \
        GET_FUNC2_MACRO(__VA_ARGS__, Func21, Func20)(__VA_ARGS__)


/*---------------------------------------------------------------------------*\
                           A templated Function3
\*---------------------------------------------------------------------------*/

class ${typeName}Function3${TemplateType}
:
    public FieldFunction3<${TemplateType}, ${typeName}Function3${TemplateType}>
{

    //- Saved dictionary
    const dictionary topDict_;
    const dictionary dict_;

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

        //- Return value as a function of scalar variable
        inline virtual ${TemplateType} value
        (
            const scalar x,
            const scalar y,
            const scalar z
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
