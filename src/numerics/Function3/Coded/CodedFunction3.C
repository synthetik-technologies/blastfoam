/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
03-12-2021 Synthetik Applied Technologies : Added Function3
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

#include "CodedFunction3.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "OSspecific.H"
#include "stringOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::Function3s::Coded<Type>::codeKeys() const
{
    return
    {
        "code",
        "codeInclude"
    };
}


template<class Type>
Foam::wordList Foam::Function3s::Coded<Type>::codeDictVars() const
{
    return {word::null, word::null};
}


template<class Type>
void Foam::Function3s::Coded<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    dynCode.setFilterVariable("typeName", codeName());

    // Set TemplateType filter variables
    dynCode.setFilterVariable("TemplateType", pTraits<Type>::typeName);

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC("codedFunction3"));

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH("codedFunction3"));

    // Make verbose if debugging
    dynCode.setFilterVariable("verbose", Foam::name(bool(debug)));

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
void Foam::Function3s::Coded<Type>::clearRedirect() const
{
    // Remove instantiation of Function3 provided by library
    redirectFunction3Ptr_.clear();
}


template<class Type>
Foam::autoPtr<Foam::Function3<Type>>
Foam::Function3s::Coded<Type>::compileNew()
{

    this->updateLibrary();

    dictionary redirectDict(dict_, codeDict());
    redirectDict.set(codeName(), codeName());

    return Function3<Type>::New(codeName(), units_, redirectDict);
}


template<class Type>
const Foam::dictionary& Foam::Function3s::Coded<Type>::expandCodeDict
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
Foam::Function3s::Coded<Type>::Coded
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    Function3<Type>(name),
    codedBase
    (
        dict.lookupOrDefault("name", name),
        expandCodeDict(dict)
    ),
    dict_(dict),
    units_(units)
{
    const fileName origCODE_TEMPLATE_DIR(getEnv("FOAM_CODE_TEMPLATES"));
    fileName tempDir(getEnv("BLAST_DIR")/"etc/codeTemplates");
    setEnv("FOAM_CODE_TEMPLATES", tempDir, true);

    redirectFunction3Ptr_ = compileNew();

    if (!origCODE_TEMPLATE_DIR.empty())
    {
        setEnv("FOAM_CODE_TEMPLATES", origCODE_TEMPLATE_DIR, true);
    }
}



template<class Type>
Foam::Function3s::Coded<Type>::Coded(const Coded<Type>& cf1)
:
    Function3<Type>(cf1),
    codedBase(cf1),
    dict_(cf1.dict_),
    units_(cf1.units_)
{
    const fileName origCODE_TEMPLATE_DIR(getEnv("FOAM_CODE_TEMPLATES"));
    fileName tempDir(getEnv("BLAST_DIR")/"etc/codeTemplates");
    setEnv("FOAM_CODE_TEMPLATES", tempDir, true);

    redirectFunction3Ptr_ = compileNew();

    if (!origCODE_TEMPLATE_DIR.empty())
    {
        setEnv("FOAM_CODE_TEMPLATES", origCODE_TEMPLATE_DIR, true);
    }
}


template<class Type>
Foam::tmp<Foam::Function3<Type>> Foam::Function3s::Coded<Type>::clone() const
{
    return tmp<Function3<Type>>(new Coded<Type>(*this));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function3s::Coded<Type>::~Coded()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function3s::Coded<Type>::value
(
    const scalarField& x,
    const scalarField& y,
    const scalarField& z
) const
{
    return
        units_.value.toStandard
        (
            redirectFunction3Ptr_->value
            (
                units_.x.toUser(x),
                units_.y.toUser(y),
                units_.z.toUser(z)
            )
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function3s::Coded<Type>::value
(
    const Field<vector>& X
) const
{
    return units_.value.toStandard
    (
        redirectFunction3Ptr_->value(units_.x.toUser(X))
    );
}

template<class Type>
void Foam::Function3s::Coded<Type>::write(Ostream& os) const
{
    writeCode(os);
}


// ************************************************************************* //
