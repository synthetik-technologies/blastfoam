/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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

#include "List.H"
#include "token.H"
#include "SLList.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::List2D<Type>::List2D(Istream& is)
:
    List<Type>(),
    m_(0),
    n_(0)
{
    operator>>(is, *this);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const List2D<Type>& L
)
{
    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII || !contiguous<Type>())
    {
        bool uniform = false;

        if (L.size() > 1 && contiguous<Type>())
        {
            uniform = true;

            forAll(L, i)
            {
                if (L.List<Type>::operator[](i) != L.List<Type>::operator[](0))
                {
                    uniform = false;
                    break;
                }
            }
        }

        if (uniform)
        {
            // Write size and start delimiter
            os << L.m_ << token::SPACE << L.n_ << token::BEGIN_BLOCK;

            // Write contents
            os << L.List<Type>::operator[](0);

            // Write end delimiter
            os << token::END_BLOCK;
        }
        else if (L.size() <= 1 || (L.n_ < 11 && contiguous<Type>()))
        {
            // Write size and start delimiter
            os  << L.m_ << token::SPACE << L.n_ << nl
                << token::BEGIN_LIST << nl;

            // Write contents
            for (label i = 0; i < L.m_; i++)
            {
                for (label j = 0; j < L.n_; j++)
                {
                    os  << L(i, j) << token::SPACE;
                }
                os  << nl;
            }

            // Write end delimiter
            os  << token::END_LIST;
        }
        else
        {
            // Write size and start delimiter
            os  << L.m_ << token::SPACE << L.n_ << nl
                << token::BEGIN_LIST << nl;

            // Write contents
            for (label i = 0; i < L.m_; i++)
            {
                for (label j = 0; j < L.n_; j++)
                {
                    os  << L(i, j) << nl;
                }
                if (i != L.m_-1) os  << nl;
            }

            // Write end delimiter
            os  << token::END_LIST;
        }
    }
    else
    {
        os << nl << L.m_ << token::SPACE << L.n_ << nl;
        if (L.size())
        {
            os.write(reinterpret_cast<const char*>(L.cdata()), L.byteSize());
        }
    }

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const List2D&)");

    return os;
}


template<class Type>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    List2D<Type>& L
)
{
    // Anull list
    L.setSize(0, 0);

    is.fatalCheck("operator>>(Istream&, List2D<Type>&)");

    token firstToken(is);

    is.fatalCheck("operator>>(Istream&, List2D<Type>&) : reading first token");

    if (firstToken.isCompound())
    {
        L.transfer
        (
            dynamicCast<token::Compound<List2D<Type>>>
            (
                firstToken.transferCompoundToken(is)
            )
        );
    }
    else if (firstToken.isLabel())
    {
        label s1 = firstToken.labelToken();
        label s2 = readLabel(is);
        label s = s1*s2;

        // Set list length to that read
        L.setSize(s1, s2);

        // Read list contents depending on data format

        if (is.format() == IOstream::ASCII || !contiguous<Type>())
        {
            // Read beginning of contents
            char delimiter = is.readBeginList("List2D");

            if (s)
            {
                if (delimiter == token::BEGIN_LIST)
                {
                    for (label i=0; i<s; i++)
                    {
                        is >> L.List<Type>::operator[](i);

                        is.fatalCheck
                        (
                            "operator>>(Istream&, List2D<Type>&) : reading entry"
                        );
                    }
                }
                else
                {
                    Type element;
                    is >> element;

                    is.fatalCheck
                    (
                        "operator>>(Istream&, List2D<Type>&) : "
                        "reading the single entry"
                    );

                    for (label i=0; i<s; i++)
                    {
                        L.List<Type>::operator[](i) = element;
                    }
                }
            }

            // Read end of contents
            is.readEndList("List2D");
        }
        else
        {
            if (s)
            {
                is.read(reinterpret_cast<char*>(L.data()), s*sizeof(Type));

                is.fatalCheck
                (
                    "operator>>(Istream&, List2D<Type>&) : reading the binary block"
                );
            }
        }
    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorInFunction(is)
                << "incorrect first token, expected '(', found "
                << firstToken.info()
                << exit(FatalIOError);
        }

        // Putback the opening bracket
        is.putBack(firstToken);

        // Now read as a singly-linked list
        SLList<List<Type>> sll(is);

        // Convert the singly-linked list to this list
        label s1 = sll.size();
        label s2 = s1 ? sll.first().size() : 0;
        L.setSize(s1, s2);
        label i = 0;
        forAllConstIter(typename SLList<List<Type>>, sll, iter)
        {
            if (s2 != iter().size())
            {
                FatalIOErrorInFunction(is)
                    << "inconsistent row sizes" << endl
                    << exit(FatalIOError);
            }
            forAll(iter(), j)
            {
                L(i, j) = iter()[j];
            }
            i++;
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int> or '(', found "
            << firstToken.info()
            << exit(FatalIOError);
    }

    return is;
}


// ************************************************************************* //
