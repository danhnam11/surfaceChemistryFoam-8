/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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

#include "porousWallHeatFlux.H"
#include "thermophysicalTransportModel.H"
#include "solidThermo.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(porousWallHeatFlux, 0);
    addToRunTimeSelectionTable(functionObject, porousWallHeatFlux, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::porousWallHeatFlux::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Wall heat-flux");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "integral");
    file() << endl;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::porousWallHeatFlux::calcWallHeatFlux
(
    const volScalarField& alpha,
    const volScalarField& he,
    const volScalarField& kappaEff,
    const volScalarField& Cp

)
{
    tmp<volScalarField> tporousWallHeatFlux
    (
        volScalarField::New
        (
            type(),
            mesh_,
            dimensionedScalar(dimMass/pow3(dimTime), 0)
        )
    );

    volScalarField::Boundary& porousWallHeatFluxBf =
        tporousWallHeatFlux.ref().boundaryFieldRef();

    const volScalarField::Boundary& heBf = he.boundaryField();
    const volScalarField::Boundary& alphaBf = alpha.boundaryField();
    const volScalarField::Boundary& kappaEffBf = kappaEff.boundaryField();
    const volScalarField::Boundary& CpBf = Cp.boundaryField();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        const label patchi = iter.key();
	if (foundObject<volScalarField>("kappaEff"))
	{
        //porousWallHeatFluxBf[patchi] = alphaBf[patchi]*heBf[patchi].snGrad();
            porousWallHeatFluxBf[patchi] = (kappaEffBf[patchi]/CpBf[patchi])*heBf[patchi].snGrad();
        }
	else
	{
	    porousWallHeatFluxBf[patchi] = alphaBf[patchi]*heBf[patchi].snGrad();
	}
    }

    if (foundObject<volScalarField>("qr"))
    {
        const volScalarField& qr = lookupObject<volScalarField>("qr");

        const volScalarField::Boundary& radHeatFluxBf = qr.boundaryField();

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            const label patchi = iter.key();

            porousWallHeatFluxBf[patchi] -= radHeatFluxBf[patchi];
        }
    }

    return tporousWallHeatFlux;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::porousWallHeatFlux::porousWallHeatFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    patchSet_()
{
    read(dict);
    resetName(typeName);
    resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::porousWallHeatFlux::~porousWallHeatFlux()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::porousWallHeatFlux::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall heat-flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::porousWallHeatFlux::execute()
{
    word name(type());

    if
    (
        foundObject<thermophysicalTransportModel>
        (
            thermophysicalTransportModel::typeName
        )
    )
    {
	
        const thermophysicalTransportModel& ttm =
            lookupObject<thermophysicalTransportModel>
            (
                thermophysicalTransportModel::typeName
            );

	if (mesh_.foundObject<volScalarField>("kappaEff"))
	{

            const volScalarField& kappaEff_(mesh_.lookupObject<volScalarField>("kappaEff"));
            return store
            (
                name,
                calcWallHeatFlux(ttm.alphaEff(), ttm.thermo().he(),kappaEff_, ttm.thermo().Cp())
            );
        }
	else
	{
	    FatalErrorInFunction
	    << "there is no kappaEff" << exit(FatalError);
	}
    }
    else if (foundObject<solidThermo>(solidThermo::dictName))
    {
        if (mesh_.foundObject<volScalarField>("kappaEff"))
	{
            const solidThermo& thermo =
                lookupObject<solidThermo>(solidThermo::dictName);
            const volScalarField& kappaEff_(mesh_.lookupObject<volScalarField>("kappaEff"));

            return store(name, calcWallHeatFlux(thermo.alpha(), thermo.he(), kappaEff_, thermo.Cp()));
	    Info << "there is kappaEff" << endl;
	}
        else
        {
            FatalErrorInFunction
            << "there is no kappaEff" << exit(FatalError);
        }

    }
    else
    {
        FatalErrorInFunction
            << "Unable to find compressible turbulence model in the "
            << "database" << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::porousWallHeatFlux::write()
{
    Log << type() << " " << name() << " write:" << nl;

    //writeLocalObjects::write();

    logFiles::write();

    const volScalarField& porousWallHeatFlux =
        obr_.lookupObject<volScalarField>(type());

    const volScalarField& kappaEff =
        obr_.lookupObject<volScalarField>("kappaEff");

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = porousWallHeatFlux.boundaryField()[patchi];

        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);
        const scalar integralHfp = gSum(magSf[patchi]*hfp);

        if (Pstream::master())
        {
            file()
                << mesh_.time().value()
                << tab << pp.name()
                << tab << minHfp
                << tab << maxHfp
                << tab << integralHfp
                << endl;
        }

        Log << "    min/max/integ(" << pp.name() << ") = "
            << minHfp << ", " << maxHfp << ", " << integralHfp << endl;
	//Log << " Jaehun Test " << " porousWallHeatFlux is : " << hfp << endl;
	//Log << " Jaehun Test " << " kappaEff : " << kappaEff << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
