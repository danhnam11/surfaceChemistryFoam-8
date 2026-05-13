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
    Foam::SurfaceReactionList

    developed by Danh Nam Nguyen
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).

Description
    List of templated reactions

SourceFiles
    SurfaceReactionList.C

\*---------------------------------------------------------------------------*/

#ifndef SurfaceReactionList_H
#define SurfaceReactionList_H

#include "SurfaceReaction.H"
#include "speciesTable.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class SurfaceReactionList Declaration
\*---------------------------------------------------------------------------*/

template<class GasThermoType, class SolidThermoType>
class SurfaceReactionList
:
    public PtrList<SurfaceReaction<GasThermoType, SolidThermoType>>
{

public:

    // Constructors

        //- Construct empty
        SurfaceReactionList()
        {}

        //- Construct from thermo list and dictionary
        SurfaceReactionList
        (
            const speciesTable& gasSpecies,
            const speciesTable& species,            
            const PtrList<GasThermoType>& gasSpeciesThermo,
            const PtrList<SolidThermoType>& speciesThermo,            
            const dictionary& dict
        );

        //- Construct from thermo list, objectRegistry and dictionary
        SurfaceReactionList
        (
            const speciesTable& gasSpecies,
            const speciesTable& species,  
            const PtrList<GasThermoType>& gasSpeciesThermo,
            const PtrList<SolidThermoType>& speciesThermo,   
            const objectRegistry& ob,
            const dictionary& dict
        );

        //- Construct copy
        SurfaceReactionList
        (
            const SurfaceReactionList<GasThermoType, SolidThermoType>& reactions
        ) = delete;


    // Public Member Functions

        //- Write
        void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SurfaceReactionList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SurfaceReactionList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
