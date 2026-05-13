/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::psiReactionCatThermo

    developed by Danh Nam Nguyen
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).

Description
    Foam::psiReactionCatThermo

SourceFiles
    psiReactionCatThermo.C

\*---------------------------------------------------------------------------*/

#ifndef psiReactionCatThermo_H
#define psiReactionCatThermo_H

#include "psiCatThermo.H"
#include "basicSpecieCatMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class psiReactionCatThermo Declaration
\*---------------------------------------------------------------------------*/

class psiReactionCatThermo
:
    public psiCatThermo
{

public:

    //- Runtime type information
    TypeName("psiReactionCatThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiReactionCatThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors

        //- Construct from mesh and phase name
        psiReactionCatThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<psiReactionCatThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~psiReactionCatThermo();


    // Member Functions

        //- Return the composition of the multi-component mixture
        virtual basicSpecieCatMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicSpecieCatMixture& composition() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
