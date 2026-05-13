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
    Foam::noSurfaceChemistrySolver

Description
    Dummy chemistry solver for 'none' option

SourceFiles
    noSurfaceChemistrySolver.H
    noSurfaceChemistrySolver.C

\*---------------------------------------------------------------------------*/

#ifndef noSurfaceChemistrySolver_H
#define noSurfaceChemistrySolver_H

#include "surfaceChemistrySolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noSurfaceChemistrySolver Declaration
\*---------------------------------------------------------------------------*/

template<class SurfaceChemistryModel>
class noSurfaceChemistrySolver
:
    public surfaceChemistrySolver<SurfaceChemistryModel>
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from thermo
        noSurfaceChemistrySolver
        (
            typename SurfaceChemistryModel::reactionThermo& thermo,
            typename SurfaceChemistryModel::reactionCatThermo& catThermo            
        );


    //- Destructor
    virtual ~noSurfaceChemistrySolver();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalar& p,
            scalar& T,
            scalarField& c,
            const label li,
            scalar& deltaT,
            scalar& subDeltaT
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "noSurfaceChemistrySolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
