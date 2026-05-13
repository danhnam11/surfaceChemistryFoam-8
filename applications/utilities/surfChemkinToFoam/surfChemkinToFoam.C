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

#include "argList.H"
#include "OFstream.H"
#include "OStringStream.H"
#include "IStringStream.H"

#include "surfChemkinReader.H"
#include "surfThermo.H"
#include "surfReactions.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    // Increase the precision of the output for JANAF coefficients
    Ostream::defaultPrecision(10);

    // Surface-phase CHEMKIN inputs
    argList::validArgs.append("CHEMKIN surface chemistry file");
    argList::validArgs.append("CHEMKIN surface thermodynamics file");
    argList::validArgs.append("surfaceThermoDict");

    // OpenFOAM outputs
    argList::validArgs.append("OpenFOAM surface chemistry file");
    argList::validArgs.append("OpenFOAM surface thermodynamics file");

    argList::addBoolOption
    (
        "newFormat",
        "read Chemkin thermo file in new format"
    );

    argList args(argc, argv);

    // Add the species thermo formatted entries
    {
        OStringStream os;
    }

    // Temporary hack to splice the specie composition data into the thermo file
    // pending complete integration into the thermodynamics structure

    for (int i=0; i<args.size(); i++)
    {
        Info << "args[" << i << "] = " << args[i] << nl;
    }

    writeSurfaceReactions(args[1], args[4]);
    //
    writeSurfaceThermo(args[1], args[2], args[5], args[3], args.optionFound("newFormat"));

    Info<< "End\n" << endl;

    Info<< nl;
    Info<< "surfChemkinToFoam execution complete." << nl;
    Info<< nl;
    Info<< "Conversion from CHEMKIN to OpenFOAM format finished successfully." << nl;
    Info<< "Output files generated:" << nl;
    Info<< "  - " << args[4] << nl;
    Info<< "  - " << args[5] << nl;
    Info<< nl;
    Info<< "Copy these files into the 'constant' directory of your simulation case before running catalystFoam." << nl;
    Info<< nl;

    return 0;
}


// ************************************************************************* //
