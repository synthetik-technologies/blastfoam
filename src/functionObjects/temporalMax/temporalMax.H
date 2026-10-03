/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           |
     \\/     M anipulation  | Synthetik Applied Technologies
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

Class
    Foam::functionObjects::temporalMax

Description
    Stores the temporal maximum value of specified field variable(s) over an
    (optionally specified) time interval of the simulation at each cell on the
    coarse mesh.

    Even if adaptive mesh refinement (AMR) dynamically refines and unrefines cells,
    each active cell in the mesh is mapped back to its coarse mesh cell ancestor,
    ensuring that the true maximum peak is preserved on the coarse mesh
    without being diluted by volume averaging during unrefinement.

    Outputs fine-grained radial peak attenuation profiles sorted
    by distance from a user-specified center.

Usage
    Example of function object specification:
    \verbatim
    radialPeakPressure
    {
        type            temporalMax;
        libs            ("libblastFunctionObjects.so");
        field           p;            // or fields ( p impulse );
        centre          (0 0 0);      // default: (0 0 0)
        timeStart       0;            // optional start time [s]
        timeEnd         1e-3;         // optional end time [s]
        writeControl    writeTime;
        outputRadial    true;         // default: true
    }
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    temporalMax.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_temporalMax_H
#define functionObjects_temporalMax_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "polyTopoChangeMap.H"
#include "HashTable.H"
#include "pointField.H"
#include "labelList.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class temporalMax Declaration
\*---------------------------------------------------------------------------*/

class temporalMax
:
    public fvMeshFunctionObject
{
protected:

    // Protected Data

        //- Names of fields to track
        wordList fieldNames_;

        //- Center coordinate for radial distance calculation
        point centre_;

        //- Start time of evaluation interval
        scalar timeStart_;

        //- End time of evaluation interval
        scalar timeEnd_;

        //- Output radial profile flag
        bool outputRadial_;

        //- Number of coarse mesh cells at initialization
        label nCoarseCells_;

        //- Centroids of the base coarse mesh cells
        pointField coarseCellCentres_;

        //- Mapping from current mesh cell index to coarse mesh cell index
        labelList cellToCoarse_;

        //- Table of accumulated maximum field values indexed by coarse cell ID
        HashTable<scalarField, word> coarseMax_;


    // Protected Member Functions

        //- Lookup or create registered volScalarField on mesh
        volScalarField& lookupOrCreate
        (
            const word& name,
            const dimensionSet& dims
        ) const;


public:

    //- Runtime type information
    TypeName("temporalMax");


    // Constructors

        //- Construct from Time and dictionary
        temporalMax
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        temporalMax(const temporalMax&) = delete;


    //- Destructor
    virtual ~temporalMax();


    // Member Functions

        //- Read the dictionary settings
        virtual bool read(const dictionary& dict);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return fieldNames_;
        }

        //- Accumulate peak values across cells at current time step
        virtual bool execute();

        //- Write out current peak values and radial profile
        virtual bool write();

        //- Handle dynamic mesh topological changes (AMR refinement/unrefinement)
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Handle mesh mapping
        virtual void mapMesh(const polyMeshMap& map);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const temporalMax&) = delete;
};

// Typedef alias for backwards compatibility
typedef temporalMax coarseMeshMax;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
