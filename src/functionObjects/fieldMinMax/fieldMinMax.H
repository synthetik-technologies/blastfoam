/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Synthetik Applied Technologies:  |  Store field maxes over all times
07-09-2021 Synthetik Applied Technologies:  |  Added minimum option
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

Class
    Foam::functionObjects::fieldMinMax

Description
    Tracks the maximum value of a given list of field names at every cell. If
    Adaptive mesh refinement is used, the maximum value may decrease due to
    volume averaging.

    Example of function object specification:
    \verbatim
    max
    {
        type            fieldMax;

        writeControl    writeTime;
        restartOnRestart    false;
        mode            component;
        minMax          min;

        fields
        (
            p
            T
        );
    }
    \endverbatim

Usage
    \table
        Property          | Description               | Required | Default
        restartOnRestart  | Restart the averaging on restart |no | no
        fields            | Name of  fields           | yes
        mode              | Mode to use               | no       | component
        minMax            | Minimum of maximum        | no       | max
    \endtable


See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    fieldMinMax.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldMinMax_H
#define functionObjects_fieldMinMax_H

#include "fvMeshFunctionObject.H"
#include "mapPolyMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{


/*---------------------------------------------------------------------------*\
                        Class fieldMax Declaration
\*---------------------------------------------------------------------------*/

class fieldMinMax
:
    public fvMeshFunctionObject
{
public:

    enum class modeType
    {
        cmpt,
        mag,
        cmptMag
    };

    enum class minMaxType
    {
        min,
        max
    };


protected:

    static const NamedEnum<modeType, 3> modeTypeNames_;
    static const NamedEnum<minMaxType, 2> minMaxTypeNames_;

    // Protected data

        //- Restart the averaging process on output
        Switch restartOnRestart_;

        //- Use magnitude or component wise
        modeType mode_;

        //- Min or max function
        minMaxType minMax_;

        //- Min max name
        word minMaxName_;

        //- Name of fields
        wordList fieldNames_;

        autoPtr<labelList> cellMap_;
        autoPtr<labelList> rCellMap_;

        //- Old field values
        HashPtrTable<regIOobject> oldFields_;

        //- Return the computed name of the min/max field
        word computedName(const word& fName) const
        {
            return
                IOobject::member(fName)
              + minMaxName_
              + IOobject::group(fName);
        }

        //- Store old fields
        template<class FieldType>
        bool map
        (
            const word& fName,
            const mapPolyMesh&
        );

        //- Update field max
        template<class FieldType>
        bool createMinMax(const word& fieldName);

        //- Create old field
        template<class FieldType>
        bool createOld(const word& fieldName);

        //- Update field max
        template<class FieldType>
        bool update(const word& fieldName);

        //- Compute the min/max using the mode
        template<class Type>
        Type compute(const Type& a, const Type& b) const;

        //- Compute the min/max using the mode
        template<class FieldType>
        void compute
        (
            FieldType& res,
            const FieldType& a,
            const FieldType& b
        ) const;

        //- Set all old fields and cell maps
        void setOldFields(const mapPolyMesh& mpm);

        //- Clear all old field and cell maps
        void clearOldFields();

public:

    //- Runtime type information
    TypeName("fieldMinMax");


    // Constructors

        //- Construct from Time and dictionary
        fieldMinMax
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        fieldMinMax(const fieldMinMax&) = delete;


    //- Destructor
    virtual ~fieldMinMax();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Update mesh and store map
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Calculate the field averages
        virtual bool execute();

        //- Write the field averages
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fieldMinMax&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldMinMaxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
