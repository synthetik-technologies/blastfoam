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
    Foam::functionObjects::temporalIntegral

Description
    Computes the temporal integral of specified field variable(s) over an
    (optionally specified) time interval of the simulation at each cell on the
    coarse mesh:
    \f[
        I(\mathbf{x}, t) = \int_{t_{\text{start}}}^t \phi^*(\mathbf{x}, \tau) \, d\tau
    \f]
    where \f$\phi^* = \max(\phi - \phi_{\text{ref}}, 0)\f$ if positiveOnly is true,
    or \f$\phi^* = (\phi - \phi_{\text{ref}})\f$ if positiveOnly is false.

    When applied to pressure with positiveOnly enabled, this computes the
    positive phase impulse of a blast wave:
    \f[
        I_+(\mathbf{x}) = \int \max(p(\mathbf{x}, t) - p_{\text{ref}}, 0) \, dt
    \f]

    Even if adaptive mesh refinement (AMR) dynamically refines and unrefines cells,
    active cells are mapped back to their base coarse cell with volume weighting,
    preserving the spatial integral without numerical diffusion from unrefinement.

    Outputs fine-grained radial integral/impulse profiles sorted
    by distance from a user-specified center.

Usage
    Example of function object specification:
    \verbatim
    radialImpulse
    {
        type            temporalIntegral;
        libs            ("libblastFunctionObjects.so");
        field           p;            // or fields ( p );
        pRef            101298;       // reference pressure / value (default: 0)
        positiveOnly    true;         // default: false (or true if pRef set)
        centre          (0 0 1.0);    // default: (0 0 0)
        timeStart       0;            // optional start time [s]
        timeEnd         5e-3;         // optional end time [s]
        writeControl    writeTime;
        outputRadial    true;         // default: true
    }
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::temporalMax

SourceFiles
    temporalIntegral.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_temporalIntegral_H
#define functionObjects_temporalIntegral_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "polyTopoChangeMap.H"
#include "polyDistributionMap.H"
#include "polyMeshMap.H"
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
                       Class temporalIntegral Declaration
\*---------------------------------------------------------------------------*/

class temporalIntegral
:
    public fvMeshFunctionObject
{
protected:

    // Protected Data

        //- Names of fields to integrate
        wordList fieldNames_;

        //- Center coordinate for radial distance calculation
        point centre_;

        //- Start time of evaluation interval
        scalar timeStart_;

        //- End time of evaluation interval
        scalar timeEnd_;

        //- Reference value to subtract: (phi - refValue)
        scalar refValue_;

        //- Integrate only positive values: max(phi - refValue, 0)
        bool positiveOnly_;

        //- Output radial profile flag
        bool outputRadial_;

        //- Number of coarse mesh cells at initialization (all processors;
        //  the coarse cells are numbered globally, so that the values follow
        //  the cells when the mesh is redistributed)
        label nCoarseCells_;

        //- Centroids of the base coarse mesh cells (all processors, master
        //  only)
        pointField coarseCellCentres_;

        //- Initial volumes of the base coarse mesh cells (all processors)
        scalarField coarseCellVolumes_;

        //- Mapping from current mesh cell index to coarse mesh cell index
        labelList cellToCoarse_;

        //- Table of accumulated integral field values indexed by coarse cell ID
        HashTable<scalarField, word> coarseIntegral_;


    // Protected Member Functions

        //- Lookup or create registered volScalarField on mesh
        volScalarField& lookupOrCreate
        (
            const word& name,
            const dimensionSet& dims
        ) const;


public:

    //- Runtime type information
    TypeName("temporalIntegral");


    // Constructors

        //- Construct from Time and dictionary
        temporalIntegral
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        temporalIntegral(const temporalIntegral&) = delete;


    //- Destructor
    virtual ~temporalIntegral();


    // Member Functions

        //- Read the dictionary settings
        virtual bool read(const dictionary& dict);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return fieldNames_;
        }

        //- Accumulate temporal integral across cells at current time step
        virtual bool execute();

        //- Write out current integral values and radial profile
        virtual bool write();

        //- Handle dynamic mesh topological changes (AMR refinement/unrefinement)
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Handle mesh mapping
        virtual void mapMesh(const polyMeshMap& map);

        //- Map the cell to coarse cell list on redistribution
        virtual void distribute(const polyDistributionMap& map);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const temporalIntegral&) = delete;
};

// Aliases
typedef temporalIntegral timeIntegral;
typedef temporalIntegral coarseMeshIntegral;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
