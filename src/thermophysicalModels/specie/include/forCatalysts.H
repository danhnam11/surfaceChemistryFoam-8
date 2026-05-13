/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#ifndef forGases_H
#define forGases_H

#include "specie.H"

#include "rhoConst.H"
#include "perfectGas.H"

#include "eConstThermo.H"
#include "hConstThermo.H"
#include "janafThermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "thermo.H"

#include "forThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forCat(Macro, Args...)                                                        \
    forThermo                                                                         \
    (                                                                                 \
        sutherlandTransport, sensibleEnthalpy, janafThermo, perfectGas, specie,       \
        Macro, Args                                                                   \
    );                                                                                \
    forThermo                                                                         \
    (                                                                                 \
        sutherlandTransport, sensibleInternalEnergy, janafThermo, perfectGas, specie, \
        Macro, Args                                                                   \
    );                                                                                \
    forThermo                                                                         \
    (                                                                                 \
        constTransport, sensibleEnthalpy, hConstThermo, rhoConst, specie,             \
        Macro, Args                                                                   \
    );                                                                                \
    forThermo                                                                         \
    (                                                                                 \
        constTransport, sensibleInternalEnergy, eConstThermo, rhoConst, specie,       \
        Macro, Args                                                                   \
    );                                                                                \

#define forCatalysts(Macro, Args...)                                                  \
    forCat(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
