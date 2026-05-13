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


#ifndef surfThermoUtils_H
#define surfThermoUtils_H

#include "dictionary.H"
#include "fileName.H"
#include "word.H"
#include "scalar.H"
#include <string>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfThermoUtils
{

std::string trimCopy(const string& s);

std::vector<std::string> splitTokens(const std::string& s);

std::string splitAlphaNumBoundaries(const std::string& s);

bool isIntegerToken(const std::string& s);

bool isScalarToken(const std::string& s);

bool isGlobalThermoTempLine(const std::string& line);

scalar stringToScalarSurface(const std::string& s);

void correctElementNameLocal(word& elementName);

scalar readSiteDensityFromSurfInp(const fileName& surfInpFile);

} // End namespace surfThermoUtils
} // End namespace Foam

#endif
// ************************************************************************* //
