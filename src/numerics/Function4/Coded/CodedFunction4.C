/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
01-11-2022 Synthetik Applied Technologies : Added Function4
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

\*---------------------------------------------------------------------------*/

#include "CodedFunction4.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "OSspecific.H"
#include "stringOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::Function4s::Coded<Type>::codeKeys() const
{
    return
    {
        "code",
        "codeInclude"
    };
}


template<class Type>
void Foam::Function4s::Coded<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    dynCode.setFilterVariable("typeName", codeName());

    // Set TemplateType filter variables
    dynCode.setFilterVariable("TemplateType", pTraits<Type>::typeName);

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC("codedFunction4"));

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH("codedFunction4"));

    // Debugging: make verbose
    if (debug)
    {
        dynCode.setFilterVariable("verbose", "true");
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
void Foam::Function4s::Coded<Type>::clearRedirect() const
{
    // Remove instantiation of Function4 provided by library
    redirectFunction4Ptr_.clear();
}


template<class Type>
Foam::autoPtr<Foam::Function4<Type>>
Foam::Function4s::Coded<Type>::compileNew()
{

    this->updateLibrary();

    dictionary redirectDict(dict_, codeDict());
    redirectDict.set(codeName(), codeName());

    return Function4<Type>::New(codeName(), redirectDict);
}


template<class Type>
const Foam::dictionary& Foam::Function4s::Coded<Type>::expandCodeDict
(
    const dictionary& cDict
) const
{
    dictionary& dict = const_cast<dictionary&>(cDict);

    verbatimString str(dict["code"]);
    stringOps::inplaceExpandEntry(str, cDict, true, true);
    dict.set(primitiveEntry("code", str));
    return cDict;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function4s::Coded<Type>::Coded
(
    const word& name,
    const dictionary& dict
)
:
    Function4<Type>(name),
    codedBase(name, expandCodeDict(dict)),
    dict_(dict)
{
    const fileName origCODE_TEMPLATE_DIR(getEnv("FOAM_CODE_TEMPLATES"));
    fileName tempDir(getEnv("BLAST_DIR")/"etc/codeTemplates");
    setEnv("FOAM_CODE_TEMPLATES", tempDir, true);

    redirectFunction4Ptr_ = compileNew();

    if (!origCODE_TEMPLATE_DIR.empty())
    {
        setEnv("FOAM_CODE_TEMPLATES", origCODE_TEMPLATE_DIR, true);
    }
}



template<class Type>
Foam::Function4s::Coded<Type>::Coded(const Coded<Type>& cf1)
:
    Function4<Type>(cf1),
    codedBase(cf1),
    dict_(cf1.dict_)
{
    const fileName origCODE_TEMPLATE_DIR(getEnv("FOAM_CODE_TEMPLATES"));
    fileName tempDir(getEnv("BLAST_DIR")/"etc/codeTemplates");
    setEnv("FOAM_CODE_TEMPLATES", tempDir, true);

    redirectFunction4Ptr_ = compileNew();

    if (!origCODE_TEMPLATE_DIR.empty())
    {
        setEnv("FOAM_CODE_TEMPLATES", origCODE_TEMPLATE_DIR, true);
    }
}


template<class Type>
Foam::tmp<Foam::Function4<Type>> Foam::Function4s::Coded<Type>::clone() const
{
    return tmp<Function4<Type>>(new Coded<Type>(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function4s::Coded<Type>::~Coded()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function4s::Coded<Type>::value
(
    const scalar t,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z
) const
{
    return redirectFunction4Ptr_->value(t, x, y, z);
}


template<class Type>
void Foam::Function4s::Coded<Type>::write(Ostream& os) const
{
    writeCode(os);
}


// ************************************************************************* //
