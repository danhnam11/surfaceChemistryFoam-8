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

#include "surfThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
namespace surfThermoUtils
{

std::string trimCopy(const string& s)
{
    std::string t(s.c_str());

    const std::string ws(" \t\r\n");
    const std::size_t first = t.find_first_not_of(ws);

    if (first == std::string::npos)
    {
        return string("");
    }

    const std::size_t last = t.find_last_not_of(ws);
    return string(t.substr(first, last - first + 1));
}

std::vector<std::string> splitTokens(const std::string& s)
{
    std::istringstream iss(s);
    std::vector<std::string> tokens;
    std::string tok;

    while (iss >> tok)
    {
        tokens.push_back(tok);
    }

    return tokens;
}

std::string splitAlphaNumBoundaries(const std::string& s)
{
    std::string out;

    for (std::size_t i = 0; i < s.size(); ++i)
    {
        const char c = s[i];
        out += c;

        if (i + 1 < s.size())
        {
            const char n = s[i + 1];

            const bool digitToAlpha =
                std::isdigit(c) && std::isalpha(n);

            const bool alphaToDigit =
                std::isalpha(c) && std::isdigit(n);

            if (digitToAlpha || alphaToDigit)
            {
                out += ' ';
            }
        }
    }

    return out;
}

bool isIntegerToken(const std::string& s)
{
    if (s.empty())
    {
        return false;
    }

    std::size_t i = 0;
    if (s[0] == '+' || s[0] == '-')
    {
        i = 1;
    }

    if (i >= s.size())
    {
        return false;
    }

    for (; i < s.size(); ++i)
    {
        if (!std::isdigit(s[i]))
        {
            return false;
        }
    }

    return true;
}

bool isScalarToken(const std::string& s)
{
    std::string t = trimCopy(s.c_str());

    if (t.empty())
    {
        return false;
    }

    bool hasDigit = false;
    bool hasExp = false;
    bool hasDot = false;

    for (std::size_t i = 0; i < t.size(); ++i)
    {
        const char c = t[i];

        if (std::isdigit(static_cast<unsigned char>(c)))
        {
            hasDigit = true;
            continue;
        }

        if ((c == '+' || c == '-') && i == 0)
        {
            continue;
        }

        if (c == '.' )
        {
            if (hasDot)
            {
                return false;
            }

            hasDot = true;
            continue;
        }

        if (c == 'e' || c == 'E' || c == 'd' || c == 'D')
        {
            if (hasExp)
            {
                return false;
            }

            hasExp = true;

            if (i + 1 < t.size() && (t[i + 1] == '+' || t[i + 1] == '-'))
            {
                ++i;
            }

            continue;
        }

        return false;
    }

    return hasDigit;
}

bool isGlobalThermoTempLine(const std::string& line)
{
    std::string t = trimCopy(line.c_str());

    if (t.empty())
    {
        return false;
    }

    std::vector<std::string> tokens = splitTokens(t);

    if (tokens.size() != 3)
    {
        return false;
    }

    return
        isScalarToken(tokens[0])
     && isScalarToken(tokens[1])
     && isScalarToken(tokens[2]);
}

scalar stringToScalarSurface(const std::string& s)
{
    std::string t = s;
    std::replace(t.begin(), t.end(), 'D', 'E');
    std::replace(t.begin(), t.end(), 'd', 'e');
    return std::atof(t.c_str());
}

void correctElementNameLocal(word& elementName)
{
    if (elementName.size() == 2)
    {
        elementName[1] = std::tolower(elementName[1]);
    }
    else if (elementName.size() == 1 && elementName[0] == 'E')
    {
        elementName = "e";
    }
}

scalar readSiteDensityFromSurfInp(const fileName& surfInpFile)
{
    IFstream is(surfInpFile);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open surface CHEMKIN file: " << surfInpFile
            << exit(FatalError);
    }

    string line;

    while (true)
    {
        is.getLine(line);

        if (!is.good())
        {
            break;
        }

        std::string s(line.c_str());

        std::size_t excl = s.find('!');
        if (excl != std::string::npos)
        {
            s.erase(excl);
        }

        if (s.empty())
        {
            continue;
        }
        if (s.find("SDEN") != std::string::npos)
        {
            std::size_t p = s.find("SDEN");
            std::size_t slash1 = s.find('/', p);
            if (slash1 == std::string::npos)
            {
                continue;
            }

            std::size_t slash2 = s.find('/', slash1 + 1);
            if (slash2 == std::string::npos)
            {
                continue;
            }

            std::string val = s.substr(slash1 + 1, slash2 - slash1 - 1);

            string foamVal(val.c_str());
            foamVal.replaceAll(" ", "");
            foamVal.replaceAll("D", "e");
            foamVal.replaceAll("d", "e");

            scalar sden = atof(foamVal.c_str());

            if (sden <= 0)
            {
                FatalErrorInFunction
                    << "Invalid SDEN value in " << surfInpFile
                    << ": " << val
                    << exit(FatalError);
            }

            Info<< "Read site density from surf.inp: " << sden << nl;
            return sden;
        }
    }

    FatalErrorInFunction
        << "No SDEN entry found in " << surfInpFile
        << exit(FatalError);

    return -1;
}

}
}
// ************************************************************************* //
