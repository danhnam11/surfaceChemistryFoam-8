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
    kineticFitChemkinToFoam

Description
    Converts CHEMKINIII thermodynamics and reaction data files into
    OpenFOAM format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "chemkinReader.H"
#include "OFstream.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "IFstream.H"
#include <sstream>

#include "fvCFD.H"
#include "fluidThermoMomentumTransportModel.H"
#include "psiReactionThermophysicalTransportModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"

#include "PPTransport.H"
#include "janafThermo.H"
#include "perfectGas.H"
#include "sensibleEnthalpy.H"
#include "DTMMultiComponentMixture.H"

#include "DPOLFT.H"
#include "DPCOEF.H"
#include "DP1VLU.H"
#include "PropertyReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

std::vector<std::string> generateTransportBlock(
    const std::string& speciesName,
    int speciesIndex,
    int nSpecies,
    const List<List<scalar>> muTC,
    const List<List<scalar>> kappaTC,
    const List<List<List<scalar>>> DijTC
)
{
    std::vector<std::string> lines;
    lines.push_back("    transport");
    lines.push_back("    {");

    std::ostringstream mu;
    mu << "        muCoeffs        (";
    for (int i = 0; i < muTC[speciesIndex].size(); ++i)
    {
        mu << muTC[speciesIndex][i];
        if (i < muTC[speciesIndex].size() - 1) mu << " ";
    }
    mu << ");";
    lines.push_back(mu.str());

    std::ostringstream kappa;
    kappa << "        kappaCoeffs     (";
    for (int i = 0; i < kappaTC[speciesIndex].size(); ++i)
    {
        kappa << kappaTC[speciesIndex][i];
        if (i < kappaTC[speciesIndex].size() - 1) kappa << " ";
    }
    kappa << ");";
    lines.push_back(kappa.str());

    lines.push_back("        DijCoeffs       (");
    for (int j = 0; j < nSpecies; ++j)
    {
        std::ostringstream dij;
        dij << "                            (";
        for (int d = 0; d < DijTC[speciesIndex][j].size(); ++d)
        {
            dij << DijTC[speciesIndex][j][d];
            if (d < DijTC[speciesIndex][j].size() - 1) dij << " ";
        }
        dij << ")";
        lines.push_back(dij.str());
    }
    lines.push_back("        );");
    lines.push_back("    }");
    return lines;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    const std::string inputFile = "constant/thermo.DTM";
    const std::string outputFile = "constant/thermo.FTM";

    std::ifstream in(inputFile);
    std::ofstream out(outputFile);

    if (!in.is_open())
    {
        std::cerr << "Failed to open " << inputFile << std::endl;
        return 1;
    }

    HashTable<label> speciesToIndex;
    forAll(speciesList, i)
    {
        speciesToIndex.insert(speciesList[i], i);
    }

    std::string line;
    std::string trimmed;
    std::string currentSpecies;
    std::string pendingSpecies;
    bool waitingForOpenBrace = false;
    bool insideElements = false;
    bool insideTransport = false;
    bool writingTransport = false;

    while (std::getline(in, line))
    {
        trimmed = line;
        trimmed.erase(0, trimmed.find_first_not_of(" \t"));
        trimmed.erase(trimmed.find_last_not_of(" \t") + 1);

        if (!waitingForOpenBrace &&
            trimmed != "" &&
            trimmed != "{" &&
            trimmed != "}" &&
            trimmed != "elements" &&
            trimmed != "specie" &&
            trimmed != "thermodynamics" &&
            trimmed != "transport" &&
            trimmed.find(' ') == std::string::npos)
        {
            pendingSpecies = trimmed;
            waitingForOpenBrace = true;
        }

        if (waitingForOpenBrace && trimmed == "{")
        {
            currentSpecies = pendingSpecies;
            waitingForOpenBrace = false;
        }

        if (trimmed == "transport")
        {
            insideTransport = true;
            writingTransport = true;
            continue;
        }

        if (insideTransport && trimmed == "}")
        {
            insideTransport = false;
            writingTransport = false;
            continue;
        }

        if (writingTransport)
        {
            continue;
        }

        if (trimmed == "elements")
        {
            insideElements = true;
        }

        if (insideElements && trimmed == "}")
        {
            insideElements = false;
            out << line << std::endl;

            if (speciesToIndex.found(currentSpecies))
            {
                label speciesIndex = speciesToIndex[currentSpecies];
                Info << "✅ inserting transport for species: " << currentSpecies << endl;

                auto transportLines = generateTransportBlock(
                    currentSpecies, speciesIndex, speciesList.size(), muTC, kappaTC, DijTC
                );

                for (const auto& l : transportLines)
                {
                    out << l << std::endl;
                }
            }
            else
            {
                Info << "⚠️  Skipping unknown species: " << currentSpecies << endl;
            }
            continue;
        }

        out << line << std::endl;
    }

    in.close();
    out.close();

    Info << "Patched thermo.DTM → thermo.FTM" << endl;
    return 0;
}

// ************************************************************************* //
