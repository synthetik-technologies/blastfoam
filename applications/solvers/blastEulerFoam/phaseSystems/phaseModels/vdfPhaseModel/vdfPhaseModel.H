/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2019-04-29 Jeff Heylmun:    Added population balance model
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
    Foam::vdfPhaseModel

Description
    Base class for a phase model that is goverened by a set of moment transport
    equations where verlocity is an internal coordinate

SourceFiles
    vdfPhaseModel.C
    newvdfPhaseModel.C
    vdfPhaseModels.C

\*---------------------------------------------------------------------------*/

#ifndef vdfPhaseModel_H
#define vdfPhaseModel_H

#include "phaseModel.H"
#include "populationBalanceModel.H"
#include "quadratureApproximations.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class vdfPhaseModel Declaration
\*---------------------------------------------------------------------------*/

class vdfPhaseModel
:
    public phaseModel
{
protected:
    // Protected data

        //- Enumeration of moment set type
        enum sizeTypeEnum
        {
            length,
            volume,
            mass
        };

        //- Enumeration of moment set type
        enum momentSetTypeEnum
        {
            volumeFractionLength,
            volumeFractionVolume,
            volumeFractionMass,
            numberDensityLength,
            numberDensityVolume,
            numberDensityMass
        };

        //- Population balance dictionary
        IOdictionary pbeDict_;

        //- Population balance model
        autoPtr<populationBalanceModel> populationBalance_;

        //- Reference to quadrature
        velocityQuadratureApproximation& quadrature_;

        //- PtrList of volume fractions
        PtrList<volScalarField> volumeFractions_;

        //- Are sigma and theta computed
        bool computeVariance_;

        //- Velocity covariance tensor
        tmp<volTensorField> sigma_;

        //- Granular temperature
        tmp<volScalarField> Theta_;

        //- Minimum diameter
        dimensionedScalar minD_;

        //- Size index
        label sizeIndex_;

        //- Moment set type
        label momentSetType_;

        //- Is the zero order moment equal to volume fraction
        bool m0VolumeFraction_;

        //- Size abscissa type
        label sizeType_;

        //- Stored labelList used to access volumeFraction index
        labelList volumeFractionMoment_;

        //- Stored labelList used to access diameter
        labelList sizeMoment_;




public:

    // Constructors
        vdfPhaseModel
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName
        );

         //- Return clone
        autoPtr<vdfPhaseModel> clone() const;


    //- Destructor
    virtual ~vdfPhaseModel();


    // Member Functions

        //- Return the number of nodes
        virtual label nNodes() const;

        //- Constant access to alpha field for nodei
        virtual tmp<volScalarField> volumeFraction(const label nodei = -1) const;

        //- Return diameter of nodei
        virtual tmp<volScalarField> d(const label nodei = -1) const;

        //- Constant access to the mean velocity
        virtual const volVectorField& U(const label nodei = -1) const;

        //- Non-const access to the mean velocity
        virtual volVectorField& U(const label nodei = -1);

        //- Deviation of velocity nodei from the mean
        virtual tmp<volVectorField> Vs(const label nodei = -1) const;

        //- Solve the underlying population balance model
        virtual void solve();

        //- Correct mean quantities
        virtual void correct();

        //- Relative transport of moments, not used in standard phase model
        virtual void relativeTransport();

        //- Average transport of moments, not used in standard phase model
        virtual void averageTransport();

        //- Solve change in size moments based on breakup and coalesence
        virtual void solveSource();

        //- Return the maximum realizable courant number
        virtual scalar realizableCo() const;

        //- Return the max courant number
        virtual scalar CoNum() const;

        //- Read phase properties dictionary
        virtual bool read()
        {
            return true;
        }

        //- Read phase properties dictionary
        virtual bool read(const bool readOK)
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
