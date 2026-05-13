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

InClass
    Foam::makeSurfaceReaction

Description
    Macros for instantiating reactions on given thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeSurfaceReaction_H
#define makeSurfaceReaction_H

#include "SurfaceReaction.H"

#include "IrreversibleSurfaceReaction.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineSurfaceReaction(nullArg, GasThermoPhysics, SolidThermoPhysics)   \
                                                                               \
    typedef SurfaceReaction<GasThermoPhysics, SolidThermoPhysics>              \
         SurfaceReaction##GasThermoPhysics##SolidThermoPhysics;                \
                                                                               \
    defineTemplateTypeNameAndDebug(                                            \
        SurfaceReaction##GasThermoPhysics##SolidThermoPhysics, 0);             \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        SurfaceReaction##GasThermoPhysics##SolidThermoPhysics,                 \
        dictionary                                                             \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        SurfaceReaction##GasThermoPhysics##SolidThermoPhysics,                 \
        objectRegistry                                                         \
    )


#define makeSurfaceReaction(SurfaceReactionType, ReactionRate,                 \
      GasThermoPhysics, SolidThermoPhysics)                                    \
                                                                               \
  typedef SurfaceReactionType<GasThermoPhysics, SolidThermoPhysics, ReactionRate>\
    SurfaceReactionType##GasThermoPhysics##SolidThermoPhysics##ReactionRate;   \
  typedef SurfaceReaction<GasThermoPhysics, SolidThermoPhysics>                \
    SurfaceReactionType##GasThermoPhysics##SolidThermoPhysics##ReactionRate##Base;\
                                                                               \
  template<>                                                                   \
  const word                                                                   \
  SurfaceReactionType##GasThermoPhysics##SolidThermoPhysics##ReactionRate::    \
  typeName                                                                     \
  (                                                                            \
      SurfaceReactionType::typeName_()                                         \
    + ReactionRate::type().capitalise()                                        \
  );                                                                           \
                                                                               \
  addToRunTimeSelectionTable                                                   \
  (                                                                            \
    SurfaceReactionType##GasThermoPhysics##SolidThermoPhysics##ReactionRate##Base,\
    SurfaceReactionType##GasThermoPhysics##SolidThermoPhysics##ReactionRate,   \
    dictionary                                                                 \
  )


#define makeIrrSurfaceReactions(ReactionRate,                                  \
          GasThermoPhysics, SolidThermoPhysics)                                \
                                                                               \
    makeSurfaceReaction                                                        \
    (                                                                          \
        IrreversibleSurfaceReaction,                                           \
        ReactionRate,                                                          \
        GasThermoPhysics,                                                      \
        SolidThermoPhysics                                                     \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //