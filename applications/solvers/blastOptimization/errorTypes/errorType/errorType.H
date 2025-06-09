/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023-2025
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
    Foam::errorType

Description

SourceFiles
    errorType.C
    errorTypeNew.C

\*---------------------------------------------------------------------------*/

#ifndef errorCriterion_H
#define errorCriterion_H

#include "runTimeSelectionTables.H"
#include "Time.H"
#include "fvMesh.H"
#include "NamedEnum.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class errorType Declaration
\*---------------------------------------------------------------------------*/

class errorType
{
public:

    // Public enumerations
    enum Reduction
    {
        MIN,
        MIN_MAG_SQR,
        MAX,
        MAX_MAG_SQR,
        AVERAGE,
        WEIGHTED_AVERAGE,
        LATEST
    };
    static const NamedEnum<Reduction, 7> timeReductionTypeNames;

    enum Target
    {
        ZERO,
        TARGET,
        LOWER_BOUND,
        UPPER_BOUND
    };
    static const NamedEnum<Target, 4> targetTypeNames;

    static label getConfigurationNo();
    static label getNConfigurations();
    static const label configuration;
    static const label nConfigurations;


protected:
// Protected data

    //- Name of the error criteria
    const word name_;

    //- Const reference to mesh
    const Time& runTime_;

    //- Name of the region
    const word regionName_;

    //- Normalise error
    bool normalise_;

    //- Current value
    scalar value_;

    // Reduction variables

        //- Time reduction method
       const Reduction timeReduction_;

        //- Integrated value
        mutable scalar intValue_;

        //- Integrated weight
        mutable scalar intWeight_;


    // Target parameters

        //- Bounding of value
        const Target targetType_;

        //- Target value
        const scalar targetValue_;

        //- Weight of error criterion
        const scalar weight_;


    // Protected member functions

        //- Return the fvMesh
        inline const fvMesh& mesh() const
        {
            return runTime_.lookupObject<fvMesh>(regionName_);
        }

        //- Reduce the current value over time
        scalar reduceValue(const scalar f);

        //- Reduce the current value over time
        static scalar reduceField
        (
            const label method,
            const scalarField& f,
            const scalarField& w
        );

        template<class T>
        static T readConfigValue
        (
            const word& name,
            const dictionary& dict
        )
        {
            const word listName = name + "s";
            if (dict.found(listName))
            {
                if (!nConfigurations)
                {
                    WarningInFunction
                        << listName << " was specified but there are no "
                        << "configurations" << endl;
                }
                List<T> vals(dict.lookup(listName));
                return vals[configuration];
            }
            else
            {
                return dict.lookup<T>(name);
            }
        }



public:

    //- Runtime type information
    TypeName("errorType");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            errorType,
            dictionary,
            (
                const Time& runTime,
                const dictionary& dict,
                const word& region
            ),
            (runTime, dict, region)
        );

    // Constructor
    errorType
    (
        const Time& runTime,
        const dictionary& dict,
        const word& region
    );

    errorType
    (
        const Time& runTime,
        const dictionary& dict,
        const word& region,
        const Reduction& timeReduction
    );

    // Selectors

        static autoPtr<errorType> New
        (
            const Time& runTime,
            const dictionary& dict,
            const word& region
        );

        //- Clone the object
        virtual autoPtr<errorType> clone() const = 0;

        // iNew class for PtrList construction
        class iNew
        {
            const Time& runTime_;
            const word region_;

        public:

            iNew(const Time& runTime, const word& region)
            :
                runTime_(runTime),
                region_(region)
            {}

            autoPtr<errorType> operator()(Istream& is) const
            {
                return errorType::New
                (
                    runTime_,
                    dictionary(word(is), dictionary::null, is),
                    region_
                );
            }
        };


    //- Destructor
    virtual ~errorType();


    // Member Functions

        //- Return the name
        inline const word& name() const
        {
            return name_;
        }

        //- Return the field name
        virtual const word& fieldName() const = 0;

        //- Return the weight
        inline scalar weight() const
        {
            return weight_;
        }

        //- Update the error
        virtual void update() = 0;

        //- Return the difference
        virtual scalar value() const
        {
            return value_;
        }

        //- Return the error
        virtual scalar error() const;

        //- Clear the errors
        virtual void clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
