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

bool parseSurfaceThermoSpeciesLine
(
    const std::string& line,
    surfaceThermoEntry& entry
)
{

//    Info<< "DEBUG: parseSurfaceThermoSpeciesLine called with line = ["
//        << line << "]" << nl;

    if (surfThermoUtils::trimCopy(line.c_str()).empty())
    {
        return false;
    }

    std::string speciesField = line.substr(0, std::min<std::size_t>(24, line.size()));
    speciesField = surfThermoUtils::trimCopy(speciesField.c_str()).c_str();

    if (speciesField.empty())
    {
        return false;
    }

    if (speciesField == "THERMO" || speciesField == "END")
    {
        return false;
    }

    if (speciesField.find("(s)") == std::string::npos)
    {
        return false;
    }

    entry.name = word(speciesField.c_str());

    std::string rest;
    if (line.size() > 24)
    {
        rest = line.substr(24);
    }
    rest = surfThermoUtils::splitAlphaNumBoundaries(rest);
    std::vector<std::string> tokens = surfThermoUtils::splitTokens(rest);

    if (tokens.size() < 5)
    {
        FatalErrorInFunction
            << "Malformed surface thermo species line for species " << entry.name << nl
            << "Line: " << line
            << exit(FatalError);
    }

    if (!tokens.empty() && surfThermoUtils::isIntegerToken(tokens.back()))
    {
        tokens.pop_back();
    }

    if (tokens.size() < 5)
    {
        FatalErrorInFunction
            << "Incomplete surface thermo species line for species " << entry.name << nl
            << "Line: " << line
            << exit(FatalError);
    }

    const std::size_t n = tokens.size();

    if (tokens[n - 4].size() != 1)
    {
        FatalErrorInFunction
            << "Could not identify phase token for species " << entry.name << nl
            << "Line: " << line
            << exit(FatalError);
    }

    entry.phase   = word(tokens[n - 4].c_str());
    entry.Tlow    = surfThermoUtils::stringToScalarSurface(tokens[n - 3]);
    entry.Thigh   = surfThermoUtils::stringToScalarSurface(tokens[n - 2]);
    entry.Tcommon = surfThermoUtils::stringToScalarSurface(tokens[n - 1]);

    DynamicList<specieElement> composition;

    for (std::size_t i = 0; i + 1 < n - 4; i += 2)
    {
        word elem(tokens[i].c_str());
        label nAtoms = std::atoi(tokens[i + 1].c_str());

        if (!elem.size() || nAtoms == 0)
        {
            continue;
        }

        surfThermoUtils::correctElementNameLocal(elem);

        specieElement se;
        se.name() = elem;
        se.nAtoms() = nAtoms;
        composition.append(se);
    }

    entry.composition = composition.shrink();

    std::string mwName(entry.name.c_str());
    std::size_t spos = mwName.find("(s)");
    if (spos != std::string::npos)
    {
        mwName.erase(spos, 3);
    }

    scalar mw = 0.0;

    for (std::size_t i = 0; i < mwName.size();)
    {
        if (!std::isupper(static_cast<unsigned char>(mwName[i])))
        {
            FatalErrorInFunction
                << "Invalid surface species name for molecular-weight parsing: "
                << entry.name
                << exit(FatalError);
        }

        std::string elemStr;
        elemStr += mwName[i++];
        if (i < mwName.size() && std::islower(static_cast<unsigned char>(mwName[i])))
        {
            elemStr += mwName[i++];
        }

        std::string numStr;
        while (i < mwName.size() && std::isdigit(static_cast<unsigned char>(mwName[i])))
        {
            numStr += mwName[i++];
        }

        label nAtoms = numStr.empty() ? 1 : atoi(numStr.c_str());

        word elem(elemStr.c_str());
        surfThermoUtils::correctElementNameLocal(elem);

        if (atomicWeights.found(elem))
        {
            mw += nAtoms*atomicWeights[elem];
        }
        else
        {
            FatalErrorInFunction
                << "Unknown element " << elem
                << " while parsing molecular weight from species name "
                << entry.name
                << exit(FatalError);
        }
    }

    entry.molWeight = mw;

    Info<< "Parsed species: " << entry.name
        << ", MW = " << entry.molWeight
        << ", composition = " << entry.composition << nl;

    return true;
}

void readSurfaceThermo
(
    const fileName& surfThermFile,
    DynamicList<surfaceThermoEntry>& thermoEntries
)
{
    IFstream is(surfThermFile);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open surface thermo file: " << surfThermFile
            << exit(FatalError);
    }

    string line;
    bool inThermoBlock = false;

    while (true)
    {
        is.getLine(line);

        if (!is.good())
        {
            break;
        }

        //Info<< "DEBUG: readSurfaceThermo line = [\"" << line << "\"]" << nl;

        std::string t = surfThermoUtils::trimCopy(line.c_str());

        if (t.empty())
        {
            continue;
        }

	if (!inThermoBlock)
	{
	    if (t == "THERMO" || t == "THERMO ALL")
	    {
	        inThermoBlock = true;
	    }
	    continue;
	}

        if (t == "END")
        {
            break;
        }
        if (surfThermoUtils::isGlobalThermoTempLine(t))
        {
            continue;
        }

        surfaceThermoEntry entry;
        if (!parseSurfaceThermoSpeciesLine(line, entry))
        {
            FatalErrorInFunction
                << "Expected surface thermo species header, but got: " << line
                << exit(FatalError);
        }

        string c1, c2, c3;
        is.getLine(c1);
        if (!is.good())
        {
            FatalErrorInFunction
                << "Unexpected EOF while reading coeff line 1 for " << entry.name
                << exit(FatalError);
        }

        is.getLine(c2);
        if (!is.good())
        {
            FatalErrorInFunction
                << "Unexpected EOF while reading coeff line 2 for " << entry.name
                << exit(FatalError);
        }

        is.getLine(c3);
        if (!is.good())
        {
            FatalErrorInFunction
                << "Unexpected EOF while reading coeff line 3 for " << entry.name
                << exit(FatalError);
        }

        //parseSurfaceThermoCoeffs(c1, c2, c3, entry);

        thermoEntries.append(entry);
    }
}


void writeSurfaceThermo
(
    const fileName& surfInpFile,
    const fileName& surfThermFile,
    const fileName& outFile,
    const fileName& surfThermoDict,
    const bool newFormat
)
{
    IFstream dictStream(surfThermoDict);

    if (!dictStream.good())
    {
        FatalErrorInFunction
            << "Cannot open thermo input dictionary: " << surfThermoDict
            << exit(FatalError);
    }

    dictionary dict(dictStream);

    IFstream is(surfThermFile);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open surface thermo file: " << surfThermFile
            << exit(FatalError);
    }

    // ---------------------------------------------------------------------
    // Current constant-property model
    // ---------------------------------------------------------------------

    // Convert cgs to SI unit
    const scalar siteDensity = 10*surfThermoUtils::readSiteDensityFromSurfInp(surfInpFile);

    dictionary defaultsDict(dict.subDict("surfaceThermoDefaults"));

    const scalar rhoDefault = readScalar(defaultsDict.lookup("rho"));
    const scalar muDefault  = readScalar(defaultsDict.lookup("mu"));
    const scalar PrDefault  = readScalar(defaultsDict.lookup("Pr"));

    const scalar CpDefault =
        //readScalar(defaultsDict.lookupOrDefault("Cp", 0.0));
        defaultsDict.lookupOrDefault<scalar>("Cp", 0.0);

    const scalar CvDefault =
        //readScalar(defaultsDict.lookupOrDefault("Cv", 0.0));
        defaultsDict.lookupOrDefault<scalar>("Cv", 0.0);

    const scalar HfDefault =
        //readScalar(defaultsDict.lookupOrDefault("Hf", 0.0));
        defaultsDict.lookupOrDefault<scalar>("Hf", 0.0);

    // Bare site enthalpy offset from your target example

    DynamicList<surfaceThermoEntry> entries;

    //readSurfaceThermo(surfThermFile, entries, readNASA);
    readSurfaceThermo(surfThermFile, entries);

    DynamicList<word> siteSpecies;
    HashSet<word> seenSiteSpecies;
    {
        IFstream surfIs(surfInpFile);

        if (!surfIs.good())
        {
            FatalErrorInFunction
                << "Cannot open surface input file: " << surfInpFile
                << exit(FatalError);
        }

        string line;
        bool inSiteBlock = false;

        while (true)
        {
            surfIs.getLine(line);

            if (!surfIs.good())
            {
                break;
            }

            std::string t(line.c_str());

            const std::size_t excl = t.find('!');
            if (excl != std::string::npos)
            {
                t.erase(excl);
            }

            t = surfThermoUtils::trimCopy(t);

            if (t.empty())
            {
                continue;
            }

            if (!inSiteBlock)
            {
                if (t.find("SITE") == 0)
                {
                    inSiteBlock = true;
                }
                continue;
            }

            if (t == "END")
            {
                break;
            }

            std::istringstream iss(t);
            std::string sp;
/*
            while (iss >> sp)
            {
                if (sp.find("(s)") != std::string::npos)
                {
                    siteSpecies.append(word(sp.c_str()));
                }
            }
*/
	    while (iss >> sp)
	    { 
	        // remove commas etc
	        sp.erase(std::remove(sp.begin(), sp.end(), ','), sp.end());
	        sp.erase(std::remove(sp.begin(), sp.end(), ';'), sp.end());
	    
	        // skip SITE metadata or occupancy tokens
	        if (sp.find('/') != std::string::npos)
	        {
	            continue;
	        }
	    
		if (sp.size() >= 3 && sp.substr(sp.size()-3) == "(s)")
	 	{
	            word wsp(sp.c_str());

    	            if (!seenSiteSpecies.found(wsp))
		    {
		        seenSiteSpecies.insert(wsp);
		        siteSpecies.append(wsp);
    		    }
		}

	    }
        }
    }

    Info<< "Parsed surface thermo species: " << entries.size() << nl;

    OFstream os(outFile);

    os  << "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        << "  =========                 |\n"
        << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n"
        << "   \\\\    /   O peration     | Website:  https://openfoam.org\n"
        << "    \\\\  /    A nd           | Version:  8\n"
        << "     \\\\/     M anipulation  |\n"
        << "-------------------------------------------------------------------------------\n"
        << "\\*---------------------------------------------------------------------------*/\n"
        << "FoamFile\n"
        << "{\n"
        << "    version     2.0;\n"
        << "    format      ascii;\n"
        << "    class       dictionary;\n"
        << "    location    \"constant\";\n"
        << "    object      thermo.catalyst.solid.const;\n"
        << "}\n"
        << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

    os << "species\n";
    os << siteSpecies.size() << "\n";
    os << "(\n";
    forAll(siteSpecies, i)
    {
        os << siteSpecies[i] << "\n";
    }
    os << ");\n\n";

    forAll(entries, i)
    {
        const surfaceThermoEntry& e = entries[i];

        os << e.name << "\n";
        os << "{\n";

        os << "    specie\n";
        os << "    {\n";
        os << "        molWeight       " << e.molWeight << ";\n";
        os << "        siteDensity     " << siteDensity << ";\n";
        os << "    }\n";

        os << "    equationOfState\n";
        os << "    {\n";
        os << "        rho             " << rhoDefault << ";\n";
        os << "    }\n";

        os << "    thermodynamics\n";
        os << "    {\n";

        os << "        Cp              " << CpDefault << ";\n";
        os << "        Hf              " << HfDefault << ";\n";
        os << "        Cv              " << CvDefault << ";\n";

        os << "    }\n";

        os << "    transport\n";
        os << "    {\n";

	os << "        mu              " << muDefault << ";\n";
	os << "        Pr              " << PrDefault << ";\n";
        os << "    }\n";

        os << "}\n\n";
    }

    os << "// ************************************************************************* //\n";
}
// ************************************************************************* //
