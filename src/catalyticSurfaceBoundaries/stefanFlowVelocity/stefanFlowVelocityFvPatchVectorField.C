/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "stefanFlowVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

#include "basicSurfaceChemistryModel.H"
#include "psiReactionThermophysicalTransportModel.H"
#include "basicThermo.H"
#include "psiReactionThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stefanFlowVelocityFvPatchVectorField::
stefanFlowVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    refValue_(p.size())
{}


Foam::stefanFlowVelocityFvPatchVectorField::
stefanFlowVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    refValue_("refValue", dict, p.size())
{
    fvPatchVectorField::operator=(refValue_*patch().nf());
}


Foam::stefanFlowVelocityFvPatchVectorField::
stefanFlowVelocityFvPatchVectorField
(
    const stefanFlowVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    refValue_(mapper(ptf.refValue_))
{
    // Note: calculate product only on ptf to avoid multiplication on
    // unset values in reconstructPar.
    fvPatchVectorField::operator=
    (
        mapper(ptf.refValue_*ptf.patch().nf())
    );
}


Foam::stefanFlowVelocityFvPatchVectorField::
stefanFlowVelocityFvPatchVectorField
(
    const stefanFlowVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    refValue_(pivpvf.refValue_)
{}


Foam::stefanFlowVelocityFvPatchVectorField::
stefanFlowVelocityFvPatchVectorField
(
    const stefanFlowVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    refValue_(pivpvf.refValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::stefanFlowVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    m(refValue_, refValue_);
}


void Foam::stefanFlowVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const stefanFlowVelocityFvPatchVectorField& tiptf =
        refCast<const stefanFlowVelocityFvPatchVectorField>(ptf);

    refValue_.rmap(tiptf.refValue_, addr);
}


void Foam::stefanFlowVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const psiReactionThermo& thermo = db().lookupObject<psiReactionThermo>
    (
        IOobject::groupName(basicThermo::dictName, IOobject::group(internalField().name()))
    );

    const basicSurfaceChemistryModel& surfaceReaction = db().lookupObject<basicSurfaceChemistryModel>
    (
        IOobject::groupName("surfaceChemistryProperties", IOobject::group(internalField().name()))
    );

    const scalarField rhow       = thermo.rho(patchi);
    const scalarField sumRgw     = surfaceReaction.momentumFluxAtCatalyticSurfaces(patchi); // [kg/(m2.s)]

    const vectorField n = patch().nf();

    // normal velocity at the wall
    vectorField Uw(patch().size(), Zero);

    forAll(Uw, facei)
    {
        const scalar rhoi = max(rhow[facei], SMALL);
        const scalar UwAbs   = -sumRgw[facei]/rhoi;

        Uw[facei] = UwAbs*n[facei];
/*
        Info<< "patch=" << patch().name()
            << " face=" << facei
            << " rho=" << rhow[facei]
            << " sumRgw=" << sumRgw[facei]
            << " UwAbs=" << UwAbs
            << " Uw=" << Uw[facei]
            << nl;
*/
    }

    fvPatchVectorField::operator=(Uw);
    fvPatchVectorField::updateCoeffs();
}


void Foam::stefanFlowVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "refValue", refValue_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        stefanFlowVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
