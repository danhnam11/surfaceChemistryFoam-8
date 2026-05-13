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

#ifndef makeSurfaceChemistrySolver_H
#define makeSurfaceChemistrySolver_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineSurfaceChemistrySolver(Model, GReactionThermo, SReactionThermo,      \
    GThermoPhysics, SThermoPhysics)                                                \
                                                                                   \
    typedef Model<GReactionThermo, SReactionThermo, GThermoPhysics, SThermoPhysics>\
        Model##GReactionThermo##SReactionThermo##GThermoPhysics##SThermoPhysics;   \
                                                                                   \
    defineTemplateTypeNameAndDebugWithName                                         \
    (                                                                              \
      Model##GReactionThermo##SReactionThermo##GThermoPhysics##SThermoPhysics,     \
      (word(                                                                       \
      Model##GReactionThermo##SReactionThermo##GThermoPhysics##SThermoPhysics::    \
        typeName_()) + "<" + GReactionThermo::typeName + ","                       \
          + SReactionThermo::typeName + "," + GThermoPhysics::typeName() + ","     \
          + SThermoPhysics::typeName() + ">"                                       \
      ).c_str(),                                                                   \
      0                                                                            \
    );

#define makeSurfaceChemistrySolver(Solver, Model, GReactionThermo, SReactionThermo,\
  GThermoPhysics, SThermoPhysics)                                                  \
                                                                                   \
  typedef Solver<Model<GReactionThermo, SReactionThermo,                           \
    GThermoPhysics, SThermoPhysics>>                                               \
    Solver##Model##GReactionThermo##SReactionThermo##GThermoPhysics##SThermoPhysics;\
                                                                                   \
  defineTemplateTypeNameAndDebugWithName                                           \
  (                                                                                \
    Solver##Model##GReactionThermo##SReactionThermo##GThermoPhysics##SThermoPhysics,\
    (                                                                              \
      word                                                                         \
      (Solver##Model##GReactionThermo##SReactionThermo##GThermoPhysics##SThermoPhysics::\
      typeName_()                                                                  \
      )                                                                            \
      + "<" +                                                                      \
      word                                                                         \
      (Model<GReactionThermo, SReactionThermo, GThermoPhysics, SThermoPhysics>::   \
      typeName_()                                                                  \
      )                                                                            \
      + "<"                                                                        \
          + word(GReactionThermo::typeName_()) + ","                               \
          + word(SReactionThermo::typeName_()) + ","                               \
          + GThermoPhysics::typeName() + ","                                       \
          + SThermoPhysics::typeName()                                             \
      + ">>"                                                                       \
    ).c_str(),                                                                     \
    0                                                                              \
  );                                                                               \
                                                                                   \
  BasicSurfaceChemistryModel<GReactionThermo, SReactionThermo>::                   \
    addthermoConstructorToTable                                                    \
    <Solver##Model##GReactionThermo##SReactionThermo##GThermoPhysics##SThermoPhysics>\
    add##Solver##Model##GReactionThermo##SReactionThermo##GThermoPhysics##SThermoPhysics##\
thermoConstructorToTable_;

#define defineSurfaceChemistrySolvers(GReactionThermo, SReactionThermo,            \
  GThermoPhysics, SThermoPhysics)                                                  \
    defineSurfaceChemistrySolver                                                   \
    (                                                                              \
        StandardSurfaceChemistryModel,                                             \
        GReactionThermo,                                                           \
        SReactionThermo,                                                           \
        GThermoPhysics,                                                            \
        SThermoPhysics                                                             \
    )

#define makeSurfaceChemistrySolvers(Solver, GReactionThermo, SReactionThermo,      \
  GThermoPhysics, SThermoPhysics)                                                  \
    makeSurfaceChemistrySolver                                                     \
    (                                                                              \
        Solver,                                                                    \
        StandardSurfaceChemistryModel,                                             \
        GReactionThermo,                                                           \
        SReactionThermo,                                                           \
        GThermoPhysics,                                                            \
        SThermoPhysics                                                             \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
