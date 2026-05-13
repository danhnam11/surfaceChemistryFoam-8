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

#include "surfReactions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static string trimCopy(const string& s)
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

static void readSiteSpecies
(
    const fileName& surfFile,
    DynamicList<word>& siteSpecies
)
{
    IFstream is(surfFile);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open surf.inp file: " << surfFile
            << exit(FatalError);
    }

    string line;
    bool inSiteBlock = false;

    while (true)
    {
        is.getLine(line);

        if (!is.good()) break;

        std::string t = trimCopy(std::string(line.c_str()));

        if (t.empty()) continue;

        // enter SITE block
        if (!inSiteBlock)
        {
            if (t.find("SITE") == 0)
            {
                inSiteBlock = true;
            }
            continue;
        }

        // exit SITE block
        if (t == "END")
        {
            break;
        }

        // tokenize species
        std::istringstream iss(t);
        std::string sp;

        while (iss >> sp)
        {
            // remove trailing commas or weird chars if needed
            //sp.erase(std::remove(sp.begin(), sp.end(), ','), sp.end());

            //siteSpecies.append(word(sp.c_str()));

            sp.erase(std::remove(sp.begin(), sp.end(), ','), sp.end());
            sp.erase(std::remove(sp.begin(), sp.end(), ';'), sp.end());
            
            if (sp.find('/') != std::string::npos)
            {
                continue;
            }
            
            if (sp.size() < 3 || sp.substr(sp.size() - 3) != "(s)")
            {
                continue;
            }
            
            siteSpecies.append(word(sp.c_str()));
        }
    }

    if (siteSpecies.empty())
    {
        FatalErrorInFunction
            << "No SITE species found in surf.inp"
            << exit(FatalError);
    }
}

static bool isNumberString(const std::string& s)
{
    if (s.empty()) return false;

    char* end = nullptr;
    std::strtod(s.c_str(), &end);

    return end != s.c_str() && *end == '\0';
}


static void removeStoichPrefix
(
    std::string& token,
    scalar& stoich
)
{
    stoich = 1.0;

    std::size_t i = 0;

    while
    (
        i < token.size()
     && (
            std::isdigit(static_cast<unsigned char>(token[i]))
         || token[i] == '.'
        )
    )
    {
        ++i;
    }

    if (i > 0)
    {
        const std::string coeff = token.substr(0, i);
        if (isNumberString(coeff))
        {
            stoich = std::atof(coeff.c_str());
            token = token.substr(i);
        }
    }
}


static void lhsReactionOrders
(
    const std::string& reaction,
    scalar& gasOrder,
    scalar& surfOrder
)
{
    gasOrder = 0.0;
    surfOrder = 0.0;

    std::size_t pos = reaction.find("=>");
    if (pos == std::string::npos)
    {
        pos = reaction.find("=");
    }

    std::string lhs =
        (pos == std::string::npos ? reaction : reaction.substr(0, pos));

    std::replace(lhs.begin(), lhs.end(), '+', ' ');

    std::istringstream iss(lhs);
    std::string sp;

    while (iss >> sp)
    {
        scalar stoich = 1.0;
        removeStoichPrefix(sp, stoich);

        if (sp.empty())
        {
            continue;
        }

        if (sp.find("(s)") != std::string::npos)
        {
            surfOrder += stoich;
        }
        else
        {
            gasOrder += stoich;
        }
    }
}


static scalar convertSurfaceArrheniusA_CGS_to_SI
(
    const scalar A_cgs,
    const scalar gasOrder,
    const scalar surfOrder
)
{
    return
        A_cgs
       *10.0
       *std::pow(1.0e-3, gasOrder)
       *std::pow(1.0e-1, surfOrder);
}


static scalar convertActivationEnergy_kJmol_to_K(const scalar Ea_kJmol)
{
    static const scalar R = 8.31446261815324; // J/mol/K
    return 1000.0*Ea_kJmol/R;
}

void writeSurfaceReactions
(
    const fileName& surfFile,
    const fileName& outFile
)
{

    DynamicList<word> siteSpecies;
    readSiteSpecies(surfFile, siteSpecies);

    IFstream is(surfFile);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open surface chemistry file: " << surfFile
            << exit(FatalError);
    }

    DynamicList<surfaceReactionEntry> rxns;

    bool inReactions = false;
    bool hasPendingReaction = false;

    surfaceReactionEntry current;

    string line;

    while (is.getLine(line))
    {
        string trimmed = trimCopy(line);

        if (trimmed.empty()) continue;
        if (trimmed[0] == '!') continue;

        if (trimmed.empty())
        {
            continue;
        }

        if (trimmed[0] == '!')
        {
            continue;
        }
        // End of REACTIONS block
        if (!inReactions)
        {
            if (trimmed.find("REACTIONS") != string::npos)
            {
                inReactions = true;
            }
            continue;
        }

        // End of REACTIONS block
        if (trimmed == "END")
        {
            if (hasPendingReaction)
            {
                rxns.append(current);
                hasPendingReaction = false;
            }
            break;
        }

        // STICK modifier line
        if (trimmed == "STICK")
        {

            if (hasPendingReaction)
            {
                current.isSKC = true;
	    
                // Sticking coefficient A is dimensionless-like.
                // Do not apply concentration-unit conversion.
                current.A = current.A_cgs;
            }
            continue;

        }

        // COV modifier line
        if (trimmed.find("COV") == 0)
        {
            if (!hasPendingReaction)
            {
                FatalErrorInFunction
                    << "COV modifier found without pending reaction: "
                    << trimmed << nl
                    << exit(FatalError);
            }

            std::string s = trimmed;
            std::replace(s.begin(), s.end(), '/', ' ');

            std::istringstream iss(s);
            std::string key, sp;
            scalar eta, mu, epsilon;

            if (!(iss >> key >> sp >> eta >> mu >> epsilon))
            {
                FatalErrorInFunction
                    << "Malformed COV modifier line: " << trimmed << nl
                    << exit(FatalError);
            }

            coverageCoeff c;
            c.species = word(sp.c_str());
            c.eta = eta;
            c.mu = mu;
	    c.epsilon = convertActivationEnergy_kJmol_to_K(epsilon);

            current.isCKC = true;
            current.cov.append(c);

            continue;
        }

        if
        (
            trimmed == "PKC"
        )
        {
            if (!hasPendingReaction)
            {
                FatalErrorInFunction
                    << "PCK modifier found without pending reaction: "
                    << trimmed << nl
                    << exit(FatalError);
            }
        
            current.isPKC = true;
            continue;
        }

        // New reaction line begins:
        // flush previous pending reaction first
        if (hasPendingReaction)
        {
            rxns.append(current);
            hasPendingReaction = false;
        }

        std::istringstream iss(trimmed.c_str());

        std::vector<std::string> tokens;
        std::string tok;

        while (iss >> tok)
        {
            tokens.push_back(tok);
        }

        if (tokens.size() < 4)
        {
            WarningInFunction
                << "Skipping malformed surface reaction line: " << trimmed << nl;
            continue;
        }

        current = surfaceReactionEntry();
        current.isSKC = false;
	current.isCKC = false;
	current.isPKC = false;
        //current.A    = atof(tokens[tokens.size()-3].c_str());
        //current.beta = atof(tokens[tokens.size()-2].c_str());
        //current.Ta   = atof(tokens[tokens.size()-1].c_str());

	const scalar A_cgs    = atof(tokens[tokens.size()-3].c_str());
	const scalar beta     = atof(tokens[tokens.size()-2].c_str());
	const scalar Ea_kJmol = atof(tokens[tokens.size()-1].c_str());

	std::string reactionRaw;

	for (std::size_t i = 0; i < tokens.size()-3; ++i)
	{
	    if (i > 0) reactionRaw += " ";
	    reactionRaw += tokens[i];
	}


	std::string reactionStd = reactionRaw;
	
	std::size_t pos = reactionStd.find("=>");
	if (pos != std::string::npos)
	{
		reactionStd.replace(pos, 2, "=");
	}
	
	reactionStd.erase
	(
		std::remove(reactionStd.begin(), reactionStd.end(), '"'),
		reactionStd.end()
	);
	
	scalar gasOrder = 0.0;
	scalar surfOrder = 0.0;
	lhsReactionOrders(reactionRaw, gasOrder, surfOrder);

	current.A_cgs     = A_cgs;
	current.gasOrder  = gasOrder;
	current.surfOrder = surfOrder;
	
	// Default: ordinary surface Arrhenius conversion
	current.A = convertSurfaceArrheniusA_CGS_to_SI
	(
		current.A_cgs,
		current.gasOrder,
		current.surfOrder
	);
	
	current.beta = beta;
	current.Ta   = convertActivationEnergy_kJmol_to_K(Ea_kJmol);
	
	
	current.reaction = reactionStd.c_str();
	
	hasPendingReaction = true;

    }

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
        << "    object      reactions.surf;\n"
        << "}\n"
        << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

    os << "reactions\n{\n";

    forAll(rxns, i)
    {
        const auto& r = rxns[i];
	    
        word reactionType;
	    
        if (r.isCKC)
        {
            reactionType =
                "irreversibleSurfaceCoverageModificationArrheniusSurface";
        }
	else if (r.isPKC)
	{
		reactionType =
			"irreversibleSurfaceTakahashiSurface";
	}
	else if (r.isSKC)
	{
		reactionType =
			"irreversibleSurfaceStickArrheniusSurface";
	}
	else
	{
		reactionType =
			"irreversibleSurfaceArrheniusSurface";
	}

        os  << "    un-named-reaction-" << i << "\n"
            << "    {\n"
            << "        type            " << reactionType << ";\n";
    
        os  << "        reaction        " << r.reaction << ";\n"
            << "        A               " << r.A << ";\n"
            << "        beta            " << r.beta << ";\n"
            << "        Ta              " << r.Ta << ";\n";
    
        if (r.isCKC)
        {
            os << "        coeffs\n";
            os << "        " << siteSpecies.size() << "\n";
            os << "        (\n";
    
            forAll(siteSpecies, si)
            {
                const word& sp = siteSpecies[si];
    
                scalar eta = 0.0;
                scalar mu = 0.0;
                scalar epsilon = 0.0;
    
                forAll(r.cov, ci)
                {
                    if (r.cov[ci].species == sp)
                    {
                        eta = r.cov[ci].eta;
                        mu = r.cov[ci].mu;
                        epsilon = r.cov[ci].epsilon;
                        break;
                    }
                }
    
                os << "            (" << sp << " ("
                   << eta << " " << mu << " " << epsilon << "))\n";
            }
    
            os << "        );\n";
        }
    
        os << "    }\n";
    }

    os << "}\n\n";
    os << "Tlow      0;\n";
    os << "Thigh     " << VGREAT << ";\n";
}
// ************************************************************************* //
