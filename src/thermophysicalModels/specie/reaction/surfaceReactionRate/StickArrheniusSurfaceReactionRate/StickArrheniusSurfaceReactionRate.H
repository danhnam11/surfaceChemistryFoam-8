/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::StickArrheniusSurfaceReactionRate

    developed by Danh Nam Nguyen
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).

Description
    Arrhenius reaction rate given by:

        k = A * T^beta * exp(-Ta/T)

SourceFiles
    StickArrheniusSurfaceReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef StickArrheniusSurfaceReactionRate_H
#define StickArrheniusSurfaceReactionRate_H

#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class StickArrheniusSurfaceReactionRate;

Ostream& operator<<(Ostream&, const StickArrheniusSurfaceReactionRate&);


/*---------------------------------------------------------------------------*\
                    Class StickArrheniusSurfaceReactionRate Declaration
\*---------------------------------------------------------------------------*/

class StickArrheniusSurfaceReactionRate
{
    // Private Data
        //scalar A_;
        //scalar beta_;
        //scalar Ta_;

        scalar a_;
        scalar b_;
        scalar c_;


public:

    // Constructors

        //- Construct from components
        inline StickArrheniusSurfaceReactionRate
        (
            const scalar a,
            const scalar b,
            const scalar c
        );

        //- Construct from dictionary
        inline StickArrheniusSurfaceReactionRate
        (
            const speciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "StickArrheniusSurface";
        }

        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        inline scalar ddT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- is sticking coefficient reaction
        inline bool isStick() const
        {
            return true;
        }

        //- is coverage dependent reaction 
        inline bool isCoverageDependence() const
        {
            return false;
        }        

        //- is takahashi
        inline bool isTakahashi() const
        {
            return false;
        }  

        //- Third-body efficiencies (beta = 1-alpha)
        //  non-empty only for third-body reactions
        //  with enhanced molecularity (alpha != 1)
        //inline const List<Tuple2<label, scalar>>& beta() const;

        //- Coverage modification coeffs
        //  non-empty only for surface coverage dependent reactions
        inline const List<Tuple2<label, scalarField>>& coverageModificationCoeffs() const;

        //- Species concentration derivative of the pressure dependent term
        //  By default this value is 1 as it multiplies the third-body term
        inline void dcidc
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalarField& dcidc
        ) const;

        //- Temperature derivative of the pressure dependent term
        //  By default this value is 0 since ddT of molecularity is approx.0
        inline scalar dcidT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const StickArrheniusSurfaceReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "StickArrheniusSurfaceReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
