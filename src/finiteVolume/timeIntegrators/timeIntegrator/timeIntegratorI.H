/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

namespace Foam
{

// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

#define defineAddFieldType(Type, Geo)                                          \
void timeIntegrator::addOldField                                               \
(                                                                              \
    const FieldName(Geo, Type)& f                                              \
) const                                                                        \
{                                                                              \
    insertOldList(OldFieldVarName(Geo, Type), f);                              \
}                                                                              \
                                                                               \
void timeIntegrator::addDeltaField                                             \
(                                                                              \
    const FieldName(Geo, Type)& f                                              \
) const                                                                        \
{                                                                              \
    insertDeltaList(DeltaFieldVarName(Geo, Type), f);                          \
}

FOR_ALL_FIELD_TYPES(defineAddFieldType, vol);
FOR_ALL_FIELD_TYPES(defineAddFieldType, surface);
FOR_ALL_FIELD_TYPES(defineAddFieldType, point);

#undef defineAddFieldType

#define defineGetFieldType(Type, Geo)   \
const PtrList<FieldName(Geo, Type)>&  \
timeIntegrator::oldFields(const FieldName(Geo, Type)& f) const  \
{   \
    return lookupOld(OldFieldVarName(Geo, Type), f);   \
}   \
    \
PtrList<FieldName(Geo, Type)>&    \
timeIntegrator::oldFields(const FieldName(Geo, Type)& f)    \
{   \
    return lookupOld(OldFieldVarName(Geo, Type), f);   \
}   \
    \
const PtrList<FieldName(Geo, Type)>&  \
timeIntegrator::deltaFields(const FieldName(Geo, Type)& f) const    \
{   \
    return lookupDelta(DeltaFieldVarName(Geo, Type), f); \
}   \
    \
PtrList<FieldName(Geo, Type)>&    \
timeIntegrator::deltaFields(const FieldName(Geo, Type)& f)  \
{   \
    return lookupDelta(DeltaFieldVarName(Geo, Type), f); \
}

FOR_ALL_FIELD_TYPES(defineGetFieldType, vol);
FOR_ALL_FIELD_TYPES(defineGetFieldType, surface);
FOR_ALL_FIELD_TYPES(defineGetFieldType, point);

#undef defineGetFieldType

} // End namespace Foam

// ************************************************************************* //
