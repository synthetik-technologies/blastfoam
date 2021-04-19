/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "particle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::particle::propertyList_ = Foam::particle::propertyList();

const std::size_t Foam::particle::sizeofPosition_
(
    offsetof(particle, facei_) - offsetof(particle, coordinates_)
);

const std::size_t Foam::particle::sizeofFields_
(
    sizeof(particle) - offsetof(particle, coordinates_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particle::particle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    mesh_(mesh),
    coordinates_(),
    celli_(-1),
    tetFacei_(-1),
    tetPti_(-1),
    facei_(-1),
    stepFraction_(0.0),
    behind_(0.0),
    nBehind_(0),
    origProc_(Pstream::myProcNo()),
    origId_(-1)
{
    if (newFormat)
    {
        if (is.format() == IOstream::ASCII)
        {
            is  >> coordinates_ >> celli_ >> tetFacei_ >> tetPti_;

            if (readFields)
            {
                is  >> facei_ >> stepFraction_ >> behind_ >> nBehind_
                    >> origProc_ >> origId_;
            }
        }
        else
        {
            if (readFields)
            {
                is.read(reinterpret_cast<char*>(&coordinates_), sizeofFields_);
            }
            else
            {
                is.read(reinterpret_cast<char*>(&coordinates_), sizeofPosition_);
            }
        }
    }
    else
    {
        compactPosition p;

        if (is.format() == IOstream::ASCII)
        {
            is >> p.position >> p.celli;

            if (readFields)
            {
                is  >> p.facei
                    >> p.stepFraction
                    >> p.tetFacei
                    >> p.tetPti
                    >> p.origProc
                    >> p.origId;
            }
        }
        else
        {
            if (readFields)
            {
                // Read whole struct
                const size_t s =
                (
                    sizeof(compactPosition)
                  - offsetof(compactPosition, position)
                );
                is.read(reinterpret_cast<char*>(&p.position), s);
            }
            else
            {
                // Read only position and cell
                const size_t s =
                (
                    offsetof(compactPosition, facei)
                  - offsetof(compactPosition, position)
                );
                is.read(reinterpret_cast<char*>(&p.position), s);
            }
        }

        if (readFields)
        {
            // Note: other position-based properties are set using locate(...)
            stepFraction_ = p.stepFraction;
            origProc_ = p.origProc;
            origId_ = p.origId;
        }

        locate
        (
            p.position,
            p.celli,
            false,
            "Particle initialised with a location outside of the mesh."
        );
    }

    // Check state of Istream
    is.check("particle::particle(const polyMesh& Istream&, bool, bool)");
}


void Foam::particle::writeCoordinates(Ostream& os) const
{
    if (os.format() == IOstream::ASCII)
    {
        os  << coordinates_
            << token::SPACE << celli_
            << token::SPACE << tetFacei_
            << token::SPACE << tetPti_;
    }
    else
    {
        os.write(reinterpret_cast<const char*>(&coordinates_), sizeofPosition_);
    }

    // Check state of Ostream
    os.check("particle::writeCoordinates(Ostream& os, bool) const");
}


void Foam::particle::writePosition(Ostream& os) const
{
    if (os.format() == IOstream::ASCII)
    {
        os  << position() << token::SPACE << celli_;
    }
    else
    {
        compactPosition p;

        const size_t s =
        (
            offsetof(compactPosition, facei)
          - offsetof(compactPosition, position)
        );

        p.position = position();
        p.celli = celli_;

        os.write(reinterpret_cast<const char*>(&p.position), s);
    }

    // Check state of Ostream
    os.check("particle::writePosition(Ostream& os, bool) const");
}


Foam::Ostream& Foam::operator<<(Ostream& os, const particle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << p.coordinates_
            << token::SPACE << p.celli_
            << token::SPACE << p.tetFacei_
            << token::SPACE << p.tetPti_
            << token::SPACE << p.facei_
            << token::SPACE << p.stepFraction_
            << token::SPACE << p.behind_
            << token::SPACE << p.nBehind_
            << token::SPACE << p.origProc_
            << token::SPACE << p.origId_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&p.coordinates_),
            particle::sizeofFields_
        );
    }

    return os;
}


// ************************************************************************* //
