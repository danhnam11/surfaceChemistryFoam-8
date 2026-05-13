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

Application
    chemkinToFoam

Description
    Converts CHEMKINIII thermodynamics and reaction data files into
    OpenFOAM format.

\*---------------------------------------------------------------------------*/


#ifndef surfThermo_H
#define surfThermo_H

#include "fileName.H"
#include <sstream>
#include "string.H"
#include <string>
#include "scalar.H"
#include "IFstream.H"
#include "OFstream.H"
#include "DynamicList.H"
#include "IOstreams.H"
#include "HashSet.H"

#include "surfThermoUtils.H"
#include "surfChemkinReader.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

struct surfaceThermoEntry
{
    word name;
    List<specieElement> composition;
    word phase;

    scalar molWeight;

    scalar Tlow;
    scalar Thigh;
    scalar Tcommon;

    FixedList<scalar, 7> highCpCoeffs;
    FixedList<scalar, 7> lowCpCoeffs;
};

struct siteCompositionInfo
{
    word siteSpeciesName;
    List<specieElement> siteComposition;
};

void writeSurfaceThermo
(
    const fileName& surfInpFile,
    const fileName& surfThermFile,
    const fileName& outFile,
    const fileName& surfThermoDict,
    const bool newFormat = false
);

bool parseSurfaceThermoSpeciesLine
(
    const std::string& line,
    surfaceThermoEntry& entry
);

void readSurfaceThermo
(
    const fileName& surfThermFile,
    DynamicList<surfaceThermoEntry>& entries
);


#endif
// ************************************************************************* //
