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

\*---------------------------------------------------------------------------*/

#include "opticallyThin.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
    defineTypeNameAndDebug(opticallyThin, 0);
    addToRadiationRunTimeSelectionTables(opticallyThin);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::opticallyThin::opticallyThin(const volScalarField& T)
:
    radiationModel(typeName, T),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    ),
    Tamb_(readScalar(coeffs_.lookup("Tamb")))
{}


Foam::radiationModels::opticallyThin::opticallyThin(const dictionary& dict, const volScalarField& T)
:
    radiationModel(typeName, dict, T),
    qr_
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), 0)
    ),
    Tamb_(readScalar(coeffs_.lookup("Tamb")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::opticallyThin::~opticallyThin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiationModels::opticallyThin::read()
{
    if (radiationModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiationModels::opticallyThin::calculate()
{
    a_ = absorptionEmission_->a();
    e_ = absorptionEmission_->e();
    E_ = absorptionEmission_->E();
}


Foam::tmp<Foam::volScalarField> Foam::radiationModels::opticallyThin::Rp() const
{
    return volScalarField::New
    (
        "Rp",
        4.0*absorptionEmission_->aCont()*constant::physicoChemical::sigma
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiationModels::opticallyThin::Ru() const
{
    tmp<volScalarField> Tenv4
    (
	new volScalarField
	(
	    IOobject
	    (
		"Tenv4",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("Tenv4", dimensionSet(0, 0, 0, 4, 0), pow(Tamb_,4.)),
            extrapolatedCalculatedFvPatchVectorField::typeName
        )
    );

    //const volScalarField::Internal& G = G_();
    const volScalarField::Internal E = absorptionEmission_->ECont()()();
    const volScalarField::Internal a = absorptionEmission_->aCont()()();
    const volScalarField::Internal T4 = Tenv4()();

    return 4.0*a*constant::physicoChemical::sigma*T4 -4.0*E;
}


// ************************************************************************* //
