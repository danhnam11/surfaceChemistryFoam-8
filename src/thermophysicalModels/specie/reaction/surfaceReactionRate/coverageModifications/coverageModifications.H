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
    Foam::coverageModifications

Description
    Third body efficiencies

SourceFiles
    coverageModificationsI.H

\*---------------------------------------------------------------------------*/

#ifndef coverageModifications_H
#define coverageModifications_H

#include "scalarList.H"
#include "speciesTable.H"
#include "scalarField.H"  //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class coverageModifications;
Ostream& operator<<(Ostream&, const coverageModifications&);


/*---------------------------------------------------------------------------*\
                  Class coverageModifications Declaration
\*---------------------------------------------------------------------------*/

class coverageModifications
//:
    //public scalarList
{
    // Private Data

        const speciesTable& species_;
        List<Tuple2<label, scalarField>> covCoeffs_;

public:

    // Constructors

        //- Construct fromdictionary
        inline coverageModifications
        (
            const speciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Calculate and return M, the concentration of the third-bodies
        //inline scalar M(const scalarList& c) const;

        //- return the coverage modification coefficients
        inline List<Tuple2<label, scalarField>> covCoeffs() const
        {
            return covCoeffs_;
        }

        //- Write to stream
        inline void write(Ostream& os) const;



    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const coverageModifications&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "coverageModificationsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
