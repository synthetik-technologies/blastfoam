/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "FunctionFieldSetType.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::FieldSetTypes::Function<Type, Patch, Mesh>::Function
(
    const fvMesh& mesh,
    const word& fieldName,
    const labelList& selectedCells,
    Istream& is,
    const bool write
)
:
    FieldSetType<Type, Patch, Mesh>(mesh, fieldName, selectedCells, is, write),
    Func1s_(0),
    funcs_(0)
{
    token equation(is);
    word eqn(equation.stringToken());

    dictionary dict(is);
    List<token::punctuationToken> ops;
    for (unsigned int i = 0; i < eqn.size(); i++)
    {
        if (eqn[i] == token::ADD)
        {
            ops.append(token::ADD);
        }
        else if (eqn[i] == token::SUBTRACT)
        {
            ops.append(token::SUBTRACT);
        }
        else if (eqn[i] == token::MULTIPLY)
        {
            ops.append(token::MULTIPLY);
        }
        else if (eqn[i] == token::DIVIDE)
        {
            ops.append(token::DIVIDE);
        }
    }

    eqn.replaceAll('+', " ");
    eqn.replaceAll('-', " ");
    eqn.replaceAll('*', " ");
    eqn.replaceAll('/', " ");
    eqn = "(" + eqn + ")";

    IStringStream eqnIs(dynamic_cast<string&>(eqn));

    tokenList coeffCmpts(eqnIs);
    labelList cmpts;
    tokenList eqnCmpts;
    wordList coeffs;
    wordList coeffNames;
    forAll(coeffCmpts, i)
    {
        const token& t(coeffCmpts[i]);
        if (t.isScalar())
        {
            eqnCmpts.append(t);
        }
        else if (t.isWord())
        {
            word w = t.wordToken();

            // create a space between variable and component
            w.replaceAll('(', ' ');

            // remove parenthesis
            w.erase(std::remove(w.begin(), w.end(), ')'), w.end());

            char xyz = w[w.size() - 1];
            w = '(' + w + ')';
            IStringStream nameIs(w);
            word name((wordList(nameIs))[0]);
            if (xyz == 'x')
            {
                cmpts.append(0);
                coeffs.append(name + "x");
                coeffNames.append(name);
            }
            else if (xyz == 'y')
            {
                cmpts.append(1);
                coeffs.append(name + "z");
                coeffNames.append(name);
            }
            else if (xyz == 'z')
            {
                cmpts.append(2);
                coeffs.append(name + "z");
                coeffNames.append(name);
            }
            eqnCmpts.append(token(coeffs.last()));
        }
    }

    //- Save underlying Function1s
    forAll(coeffNames, i)
    {
        if (!Func1s_.found(coeffNames[i]))
        {
            Func1s_.resize(Func1s_.size() + 1);
            Func1s_.set
            (
                i,
                coeffNames[i],
                Function1<Type>::New(coeffNames[i], dict).ptr()
            );
        }
    }

    //- Set functions with components
    forAll(coeffs, i)
    {
        if (!funcs_.found(coeffs[i]))
        {
            funcs_.resize(funcs_.size() + 1);
            funcs_.set
            (
                i,
                coeffs[i],
                new directionalFunction1<Type>(Func1s_[coeffNames[i]], cmpts[i])
            );
        }
    }

    if (ops.size() == eqnCmpts.size())
    {
        forAll(ops, i)
        {
            equation_.append(ops[i]);
            equation_.append(eqnCmpts[i]);
        }
    }
    else
    {
        equation_.append(eqnCmpts[0]);
        forAll(ops, i)
        {
            equation_.append(ops[i]);
            equation_.append(eqnCmpts[i+1]);
        }
    }

    if (this->good_)
    {
        setField();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::FieldSetTypes::Function<Type, Patch, Mesh>::~Function()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
void Foam::FieldSetTypes::Function<Type, Patch, Mesh>::setField()
{

    Field<Type> values(this->selectedCells_.size(), Zero);
    label nParts = 1;
    forAll(equation_, i)
    {
        if (equation_[i] == token::ADD || equation_[i] == token::SUBTRACT)
        {
            nParts++;
        }
    }

    List<Field<Type>> parts
    (
        nParts,
        Field<Type>(this->selectedCells_.size(), Zero)
    );
    Type one = pTraits<Type>::one;
    const vectorField& xyz = this->mesh_.C().primitiveField();
    label parti = 0;
    label i = 0;

    if (!equation_[i].isPunctuation())
    {
        if (equation_[i].isScalar())
        {
            parts[parti] = pTraits<Type>::one*equation_[i++].scalarToken();
        }
        else
        {
            parts[parti] = funcs_[equation_[i++].wordToken()].value(xyz);
        }
    }

    while (i < equation_.size())
    {
        token::punctuationToken op = equation_[i++].pToken();
        Field<Type> val(values.size(), one);
        if (equation_[i].isScalar())
        {
            val *= equation_[i].scalarToken();
        }
        else if (equation_[i].isWord())
        {
            val = funcs_[equation_[i].wordToken()].value(xyz);
        }

        if (op == token::ADD)
        {
            parts[++parti] = val;
        }
        else if (op == token::SUBTRACT)
        {
            parts[++parti] = -val;
        }
        else if (op == token::MULTIPLY)
        {
            parts[parti] = cmptMultiply(parts[parti], val);
        }
        else if (op == token::DIVIDE)
        {
            parts[parti] = cmptDivide(parts[parti], val);
        }
        i++;
    }
    forAll(parts, parti)
    {
        values += parts[parti];
    }

    forAll(this->selectedCells_, i)
    {
        (*this->fieldPtr_)[this->selectedCells_[i]] = values[i];
    }

    FieldSetType<Type, Patch, Mesh>::setField();
}

