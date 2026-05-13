/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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

#include "basicSurfaceChemistryModel.H"
#include "basicThermo.H"
#include "basicCatThermo.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class SurfaceChemistryModel>
Foam::autoPtr<SurfaceChemistryModel> Foam::basicSurfaceChemistryModel::New
(
    typename SurfaceChemistryModel::reactionThermo& thermo,
    typename SurfaceChemistryModel::reactionCatThermo& catThermo    
)
{
    IOdictionary chemistryDict
    (
        IOobject
        (
            thermo.phasePropertyName("surfaceChemistryProperties"),
            thermo.db().time().constant(),
            thermo.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    if (!chemistryDict.isDict("surfaceChemistryType"))
    {
        FatalErrorInFunction
            << "Template parameter based chemistry solver selection is no "
            << "longer supported. Please create a chemistryType dictionary"
            << "instead." << endl << endl << "For example, the entry:" << endl
            << "    chemistrySolver ode<StandardSurfaceChemistryModel<"
            << "rhoSurfaceChemistryModel,sutherland<specie<janaf<perfectGas>,"
            << "sensibleInternalEnergy>>>>" << endl << endl << "becomes:"
            << endl << "    chemistryType" << endl << "    {" << endl
            << "        solver ode;" << endl << "        method standard;"
            << endl << "    }" << exit(FatalError);
    }

    const dictionary& chemistryTypeDict =
        chemistryDict.subDict("surfaceChemistryType");

    const word& solverName
    (
        chemistryTypeDict.found("solver")
      ? chemistryTypeDict.lookup("solver")
      : chemistryTypeDict.found("chemistrySolver")
      ? chemistryTypeDict.lookup("chemistrySolver")
      : chemistryTypeDict.lookup("solver") // error if neither entry is found
    );

    const word& methodName
    (
        chemistryTypeDict.lookupOrDefault<word>
        (
            "method",
            chemistryTypeDict.lookupOrDefault<bool>("TDAC", false)
          ? "TDAC"
          : "standardSCHEM"
        )
    );

    dictionary chemistryTypeDictNew;
    chemistryTypeDictNew.add("solver", solverName);
    chemistryTypeDictNew.add("method", methodName);

    Info<< "Selecting surface chemistry solver " << chemistryTypeDictNew << endl;

    typedef typename SurfaceChemistryModel::thermoConstructorTable cstrTableType;
    cstrTableType* cstrTable = SurfaceChemistryModel::thermoConstructorTablePtr_;

    const word chemSolverCompThermoName =
        solverName + '<' + methodName + '<'
      + SurfaceChemistryModel::reactionThermo::typeName + ','
      + SurfaceChemistryModel::reactionCatThermo::typeName + ','
      + thermo.thermoName() + ','             
      + catThermo.thermoName() + ">>";

    typename cstrTableType::iterator cstrIter =
        cstrTable->find(chemSolverCompThermoName);

    if (cstrIter == cstrTable->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName_() << " type " << solverName << endl
            << endl;

        const wordList names(cstrTable->toc());

        // this is for gas-phase thermo physics
        wordList thisCmpts;
        thisCmpts.append(word::null);
        thisCmpts.append(word::null);
        thisCmpts.append(SurfaceChemistryModel::reactionThermo::typeName);
        thisCmpts.append(SurfaceChemistryModel::reactionCatThermo::typeName);            
        thisCmpts.append(basicThermo::splitThermoName(thermo.thermoName(), 5));
        thisCmpts.append(basicThermo::splitThermoName(catThermo.thermoName(), 5));

        List<wordList> validNames;
        validNames.append(wordList(2, word::null));
        validNames[0][0] = "solver";
        validNames[0][1] = "method";
        forAll(names, i)
        {
            const wordList cmpts(basicThermo::splitThermoName(names[i], 14));

            bool isValid = true;
            for (label i = 2; i < cmpts.size() && isValid; ++ i)
            {
                isValid = isValid && cmpts[i] == thisCmpts[i];
            }

            if (isValid)
            {
                validNames.append(SubList<word>(cmpts, 2));
            }
        }

        FatalErrorInFunction
            << "All " << validNames[0][0] << '/' << validNames[0][1]
            << " combinations for this thermodynamic model are:"
            << endl << endl;
        printTable(validNames, FatalErrorInFunction);

        FatalErrorInFunction << endl;

        List<wordList> validCmpts;
        validCmpts.append(wordList(14, word::null));
        validCmpts[0][0] = "solver";
        validCmpts[0][1] = "method";
        validCmpts[0][2] = "gReactionThermo";
        validCmpts[0][3] = "sReactionThermo";        
        validCmpts[0][4] = "transport";
        validCmpts[0][5] = "thermo";
        validCmpts[0][6] = "equationOfState";
        validCmpts[0][7] = "specie";
        validCmpts[0][8] = "energy";
        validCmpts[0][9] = "solidTransport";
        validCmpts[0][10] = "solidThermo";
        validCmpts[0][11] = "solidEquationOfState";
        validCmpts[0][12] = "solidSpecie";
        validCmpts[0][13] = "solidEnergy";

        forAll(names, i)
        {
            validCmpts.append(basicThermo::splitThermoName(names[i], 14));
        }

        FatalErrorInFunction
            << "All " << validCmpts[0][0] << '/' << validCmpts[0][1] << '/'
            << validCmpts[0][2] << '/' << validCmpts[0][3] 
            << "/thermoPhysics combinations are:" << endl << endl;
        printTable(validCmpts, FatalErrorInFunction);

        FatalErrorInFunction << exit(FatalError);
    }

    return autoPtr<SurfaceChemistryModel>(cstrIter()(thermo, catThermo));
}

// ************************************************************************* //
