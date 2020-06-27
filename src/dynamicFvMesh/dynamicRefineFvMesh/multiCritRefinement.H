/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::dynamicRefineFvMesh

Description
    Implementation of multi-criterion refinement options for dynamicRefineFvMesh

    Multiples criteria can be selected:
    fields      (e.g. alpha, speciesConcentration)
    gradients   (e.g. grad(U))
    curls       (e.g. )
    interface   (e.g. alpha 0.5 iso-surface)

    To be used with 
        dynamicRefineFvMeshCoeffs.field multiCritRefinementField
    and 
        multiCritRefinementControls
    Has no effect if the multiCritRefinementControls subdict doesnt exist.

SourceFiles
    multiCritRefinement.C

Authors
    Daniel Deising <deising@mma.tu-darmstadt.de>
    Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de>
    All rights reserved.

Description
    multi-criterion refinement as plugin for dynamicRefineFvMesh

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Daniel Rettenmaier (main developer).
    
    Method Development and Intellectual Property :
        Daniel Deising <deising@mma.tu-darmstadt.de>
        Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de>
        Holger Marschall <marschall@csi.tu-darmstadt.de>
        Dieter Bothe <bothe@csi.tu-darmstadt.de>
        Cameron Tropea <ctropea@sla.tu-darmstadt.de>

        School of Aeronautics and Astronautics
        Purdue University
        Mathematical Modeling and Analysis
        Institute for Fluid Mechanics and Aerodynamics
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.


\*---------------------------------------------------------------------------*/

#ifndef multiCritRefinement_H
#define multiCritRefinement_H

//#include "dynamicFvMesh.H"
//#include "PackedBoolList.H"
#include "dynamicFvMesh.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "syncTools.H"
#include "pointFields.H"
#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "pointMesh.H"
#include "Switch.H"
#include "topoSetSource.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dynamicRefineFvMesh Declaration
\*---------------------------------------------------------------------------*/

class multiCritRefinement
{

private:

        //- const reference to current refinement Level
        const labelList& cellLevel_;

        //- const reference to mesh
        const fvMesh& mesh_;

        //- difference of target and current refinement level > 0 refineme
        //                                                    < 0 unrefine
        volScalarField* multiCritRefinementFieldPtr_;

        //- target refinement level combining all criteria
        volScalarField* targetLevelPtr_;

        labelList targetLevel_;

        //- current refinement level
        volScalarField* isLevelPtr_;

        //- list of all fields and their refinement criteria:
        //   minValue <= field <= maxValue
        //   refinementLevel, nAddLayers                                                    
        HashTable< dictionary > fields_;
        
        //- list of scalar fields of which the |gradient| serves as refinement criteria
        //  minValue <= mag(grad(field)) <= maxValue
        HashTable< dictionary > gradFields_;

        //- list of fields of which its curl is used as refinement criteria
        //  minValue <= mag(curl(field)) <= maxValue
        HashTable< dictionary > curlFields_;

        //- inner and outer regions of an interface ( 0.1 < grad(alpha) ) are set
        //  to the global maxRefineLevel (in most cases you want exactly that)
        HashTable< dictionary > interface_;

        //- regions can be set to be quasi-static refined
        PtrList<entry> refinedRegions_;

        //- switch all criteria on or of
        Switch enableMultiCritRefinementControl_;

        //- the number of buffer layers between three different refinement levels
        //  preferrably an even number
        label nBufferLayers_;

        //- reads and loads all "multiCritRefinementControl" - information
        bool readMultiCritRefinementDict();

        //- calculates target level for the criteria "field", "gradient" and "curl"
        //  evaluates dict entries:
        //      fieldName, minValue, maxValue, refineLevel, nAddLayers
        void applyCritEntries(word critType, dictionary critDict, word critName);

        //- Disallow default bitwise assignment
        void operator=(const multiCritRefinement&);


public:

    // Constructors

        //- Constructor
            explicit multiCritRefinement(
                    const labelList& cellLevel, 
                    const fvMesh& mesh
            );


        //- Destructor
        virtual ~multiCritRefinement();


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        multiCritRefinement(const multiCritRefinement&);

        //- calculates the multiCritRefinementField 
        void updateRefinementField();



    // Writing

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
