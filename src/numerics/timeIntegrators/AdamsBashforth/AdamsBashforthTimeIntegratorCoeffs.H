/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2023
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

Class
    Foam::timeIntegrators::AdamsBashforth

SourceFiles
    AdamsBashforthTimeIntegratorCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef AdamsBashforthTimeIntegratorCoeffs_H
#define AdamsBashforthTimeIntegratorCoeffs_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegratorCoeffs.H"
#include "IntegratorsFwd.H"

namespace Foam
{
namespace timeIntegrators
{

/*---------------------------------------------------------------------------*\
                           Class AdamsBashforth Declaration
\*---------------------------------------------------------------------------*/

class AdamsBashforth
:
    public timeIntegratorCoeffs
{
public:
    // Equation for determining coefficients of arbitrary accuracy and spacing
    class LagrangePolyCoeffEqn
    :
        public ScalarEquation
    {
        mutable label index_;
        const List<scalar>& times_;

    public:

        LagrangePolyCoeffEqn(const List<scalar>& times)
        :
            index_(-1),
            times_(times)
        {}

        void setIndex(const label index) const
        {
            index_ = index;
        }
        label nDerivatives() const
        {
            return 0;
        }
        scalar fx(const scalar t, const label li) const
        {
            scalar val = 1.0;
            for (label i = 1; i < times_.size(); i++)
            {
                if (i != index_)
                {
                    val *= (t - times_[i])/(times_[index_] - times_[i]);
                }
            }
            return val;
        }
    };

private:
    //- Order of accuracy
    label order_;

    //- Current order of accuracy
    label currOrder_;

    // Saved times
    List<scalar> times_;

    //- Are the times uniformly spaced
    bool uniform_;

    //- Starting index
    mutable label startIndex_;

    //- Equation for coefficient integration
    LagrangePolyCoeffEqn eqn_;

    //- Integrator for computing weights
    autoPtr<scalarIntegrator> integrator_;


public:

    //- Runtime type information
    TypeName("AdamsBashforth");

    // Constructor
    AdamsBashforth(Istream& is);


    //- Destructor
    virtual ~AdamsBashforth();


    // Member Functions

        // Coefficient functions

            //- Set coefficients
            virtual void set
            (
                List<List<scalar>>& as,
                List<List<scalar>>& bs,
                const label index
            ) const;

            //- Update coefficients
            virtual bool update(const label index, const scalar t);

            //- Return the index of saved old fields
            virtual List<label> oldIs(const List<List<scalar>>& as) const;

            //- Return the index of saved delta fields
            virtual List<label> deltaIs(const List<List<scalar>>& bs) const;

            // Set the initial and end time coefficients (all based on current time)
            void setTimeFactors
            (
                const List<List<scalar>>& as,
                const List<List<scalar>>& bs,
                List<scalar>& f0,
                List<scalar>& f
            ) const;

            //- Return max Courant number
            virtual scalar maxCo() const
            {
                return 1.0/pow(2.0, currOrder_ - 1);
            }


        // Mapping/save functions

            //- Are there any old times that need to be saved
            virtual bool save() const
            {
                return true;
            }

            //- Return the map for saving old values
            virtual List<label> deltaSaveMap() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace timeIntegrators
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
