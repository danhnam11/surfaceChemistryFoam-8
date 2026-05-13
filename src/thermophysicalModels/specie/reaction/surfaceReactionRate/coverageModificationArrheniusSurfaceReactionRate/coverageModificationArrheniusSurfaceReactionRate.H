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
    Foam::coverageModificationArrheniusSurfaceReactionRate

    developed by Danh Nam Nguyen 
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).


Description
    Arrhenius reaction rate enhanced by third-body interaction.

SourceFiles
    coverageModificationArrheniusSurfaceReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef coverageModificationArrheniusSurfaceReactionRate_H
#define coverageModificationArrheniusSurfaceReactionRate_H

#include "ArrheniusSurfaceReactionRate.H"
#include "coverageModifications.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class coverageModificationArrheniusSurfaceReactionRate;

Ostream& operator<<(Ostream&, const coverageModificationArrheniusSurfaceReactionRate&);


/*---------------------------------------------------------------------------*\
              Class coverageModificationArrheniusSurfaceReactionRate Declaration
\*---------------------------------------------------------------------------*/

class coverageModificationArrheniusSurfaceReactionRate
:
    ArrheniusSurfaceReactionRate
{
    // Private Data

        coverageModifications coverageModifications_;
        List<Tuple2<label, scalarField>> coverageModificationCoeffs_;


public:

    // Constructors

        //- Construct from components
        inline coverageModificationArrheniusSurfaceReactionRate
        (
            const scalar A,
            const scalar beta,
            const scalar Ta,
            const coverageModifications& covM
        );

        //- Construct from dictionary
        inline coverageModificationArrheniusSurfaceReactionRate
        (
            const speciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "coverageModificationArrheniusSurface";
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
            return false;
        }

        //- is coverage dependent reaction 
        inline bool isCoverageDependence() const
        {
            return true;
        }

        //- is Takahashi
        inline bool isTakahashi() const
        {
            return false;
        }

        //- Third-body efficiencies (beta = 1-alpha)
        //  non-empty only for third-body reactions
        //  with enhanced molecularity (alpha != 1)
        /*
        inline const List<Tuple2<label, scalar>>& beta() const
        {
            return beta_;
        }
        */

        //- Coverage modification coeffs
        //  non-empty only for surface coverage dependent reactions
        inline const List<Tuple2<label, scalarField>>& coverageModificationCoeffs() const
        {
            return coverageModificationCoeffs_;
        }


        //- Species concentration derivative of the pressure dependent term
        inline void dcidc
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalarField& dcidc
        ) const;

        //- Temperature derivative of the pressure dependent term
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
            const coverageModificationArrheniusSurfaceReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "coverageModificationArrheniusSurfaceReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
