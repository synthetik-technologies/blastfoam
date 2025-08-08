/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
12-04-2022 Synthetik Applied Technologies : Added equation functionality
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "CodedCoefficientEquation.H"
#include "stringOps.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::CodedCoefficientEquation<Type>::codeKeys() const
{
    return
    {
        "fx_code",
        "dfdx_code",
        "coeffJ_code",
        "codeInclude"
    };
}


template<class Type>
Foam::wordList Foam::CodedCoefficientEquation<Type>::codeDictVars() const
{
    return
    {
        word::null,
        word::null,
        word::null,
        word::null
    };
}


template<class Type>
void Foam::CodedCoefficientEquation<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    dynCode.setFilterVariable("typeName", codeName());

    dynCode.setFilterVariable("nDerivatives", Foam::name(nDerivatives_));

    // Set TemplateType filter variables
    dynCode.setFilterVariable("TemplateType", pTraits<Type>::typeName);

    // Make verbose if debugging
    dynCode.setFilterVariable("verbose", Foam::name(bool(debug)));

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC("CodedCoefficientEquation"));

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH("CodedCoefficientEquation"));

    // Debugging: make verbose
    if (debug)
    {
        Info<<"compile " << codeName() << " sha1: "
            << context.sha1() << endl;
    }

    // Define Make/options
    dynCode.setMakeOptions
    (
        "EXE_INC = -g \\\n"
        + word("    -I$(BLAST_DIR)/src/numerics/lnInclude \\\n")
        + context.options()
        + "\n\nLIB_LIBS = \\\n"
        + "    -lOpenFOAM \\\n"
        + "    -L$(BLAST_LIBBIN) \\\n"
        + "    -lblastNumerics \\\n"
        + context.libs()
    );
}


template<class Type>
void Foam::CodedCoefficientEquation<Type>::clearRedirect() const
{
    // Remove instantiation of Equation provided by library
    redirectEquationPtr_.clear();
}


template<class Type>
Foam::autoPtr<Foam::UnivariateCoefficientEquation<Type>>
Foam::CodedCoefficientEquation<Type>::compileNew()
{
    this->updateLibrary();

    return autoPtr<UnivariateCoefficientEquation<Type>>
    (
        dynamic_cast<UnivariateCoefficientEquation<Type>*>
        (
            UnivariateEquation<Type>::New
            (
                codeName(),
                codeDict()
            ).ptr()
        )
    );
}


template<class Type>
const Foam::dictionary&
Foam::CodedCoefficientEquation<Type>::expandCodeDict
(
    const dictionary& cDict
) const
{
    const wordList codes
    (
        {
            "fx_code",
            "dfdx_code",
            "coeffJ_code",
        }
    );
    dictionary& dict = const_cast<dictionary&>(cDict);
    forAll(codes, i)
    {
        verbatimString str(dict[codes[i]]);
        stringOps::inplaceExpandEntry(str, cDict, true, true);
        dict.set(primitiveEntry(codes[i], str));
    }
    return cDict;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CodedCoefficientEquation<Type>::CodedCoefficientEquation
(
    const dictionary& dict
)
:
    UnivariateCoefficientEquation<Type>
    (
        dict,
        dict.lookup<label>("nCoeffs")
    ),
    codedBase("test", expandCodeDict(dict)),
    nDerivatives_(dict.lookup<label>("nDerivatives"))
{
    const fileName origCODE_TEMPLATE_DIR(getEnv("FOAM_CODE_TEMPLATES"));
    fileName tempDir(getEnv("BLAST_DIR")/"etc/codeTemplates");
    setEnv("FOAM_CODE_TEMPLATES", tempDir, true);

    redirectEquationPtr_ = compileNew();

    if (!origCODE_TEMPLATE_DIR.empty())
    {
        setEnv("FOAM_CODE_TEMPLATES", origCODE_TEMPLATE_DIR, true);
    }
}

template<class Type>
Foam::CodedCoefficientEquation<Type>::CodedCoefficientEquation
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    CodedCoefficientEquation<Type>(dict)
{
    this->setObr(obr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::CodedCoefficientEquation<Type>::~CodedCoefficientEquation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
