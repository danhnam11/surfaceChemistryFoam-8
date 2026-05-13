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

#ifndef forSCHEM_H
#define forSCHEM_H

#include "specie.H"

#include "perfectGas.H"
#include "rhoConst.H"

#include "janafThermo.H"
#include "eConstThermo.H"
#include "hConstThermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"

#include "sutherlandTransport.H"
#include "DTMTransport.H"
#include "FTMTransport.H"
#include "PPTransport.H"
#include "constTransport.H"

#include "thermo.H"

#include "forThermo.H"
#include "forThermoSCHEM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define forSuGasEqns(gMu, He, Cp, Macro, Args...)                          \
     forThermoSCHEM(gMu, He, Cp, perfectGas, specie,                       \
    sutherlandTransport, He, Cp, perfectGas, specie,                       \
	Macro, Args) 

#define forSuThermos(gMu, Macro, Args...)                                  \
    forSuGasEqns(gMu, sensibleEnthalpy, janafThermo, Macro, Args);         \
    forSuGasEqns(gMu, sensibleInternalEnergy, janafThermo, Macro, Args) 

#define forConstGasEqns(gMu, He, gCp, sCp, Macro, Args...)                 \
     forThermoSCHEM(gMu, He, gCp, perfectGas, specie,                      \
         constTransport, He, sCp, rhoConst, specie,                        \
	Macro, Args) 

#define forConstThermos(gMu, Macro, Args...)                               \
    forConstGasEqns(gMu, sensibleEnthalpy, janafThermo, hConstThermo,      \
	Macro, Args);                                                      \
    forConstGasEqns(gMu, sensibleInternalEnergy, janafThermo, eConstThermo,\
	Macro, Args) 

#define forSCHEMtransports(Macro, Args...)                                 \
    forSuThermos(sutherlandTransport, Macro, Args);                        \
    forSuThermos(DTMTransport, Macro, Args);                               \
    forSuThermos(FTMTransport, Macro, Args);                               \
    forSuThermos(PPTransport, Macro, Args);                                \
    forConstThermos(sutherlandTransport, Macro, Args);                     \
    forConstThermos(DTMTransport, Macro, Args);                            \
    forConstThermos(FTMTransport, Macro, Args);                            \
    forConstThermos(PPTransport, Macro, Args)

#define forSCHEM(Macro, Args...)                                           \
    forSCHEMtransports(Macro, Args)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
