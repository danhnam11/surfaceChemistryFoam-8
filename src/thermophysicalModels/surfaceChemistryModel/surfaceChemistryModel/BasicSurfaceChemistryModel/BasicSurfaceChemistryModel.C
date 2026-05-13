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

\*---------------------------------------------------------------------------*/

#include "BasicSurfaceChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class GasReactionThermo, class SolidReactionThermo>
Foam::BasicSurfaceChemistryModel<GasReactionThermo, SolidReactionThermo>::
BasicSurfaceChemistryModel
(
    GasReactionThermo& thermo,
    SolidReactionThermo& catThermo    
)
:
    basicSurfaceChemistryModel(thermo, catThermo),
    thermo_(thermo),
    catThermo_(catThermo)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class GasReactionThermo, class SolidReactionThermo>
Foam::autoPtr
<
    Foam::BasicSurfaceChemistryModel<GasReactionThermo, SolidReactionThermo>
>
Foam::BasicSurfaceChemistryModel<GasReactionThermo, SolidReactionThermo>::
New
(
    GasReactionThermo& thermo,
    SolidReactionThermo& catThermo    
)
{
    return basicSurfaceChemistryModel::
    New
    <
        BasicSurfaceChemistryModel<GasReactionThermo, SolidReactionThermo>
    >
    (
        thermo,
        catThermo
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class GasReactionThermo, class SolidReactionThermo>
Foam::BasicSurfaceChemistryModel<GasReactionThermo, SolidReactionThermo>::
~BasicSurfaceChemistryModel()
{}


// ************************************************************************* //
