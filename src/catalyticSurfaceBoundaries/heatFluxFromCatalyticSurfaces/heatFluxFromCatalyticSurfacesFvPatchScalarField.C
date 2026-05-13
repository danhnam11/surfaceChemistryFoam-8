/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------


\*---------------------------------------------------------------------------*/

#include "heatFluxFromCatalyticSurfacesFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicSurfaceChemistryModel.H"
#include "psiReactionThermophysicalTransportModel.H"
#include "basicThermo.H"
#include "psiReactionThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatFluxFromCatalyticSurfacesFvPatchScalarField::
heatFluxFromCatalyticSurfacesFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


Foam::heatFluxFromCatalyticSurfacesFvPatchScalarField::
heatFluxFromCatalyticSurfacesFvPatchScalarField
(
    const heatFluxFromCatalyticSurfacesFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::heatFluxFromCatalyticSurfacesFvPatchScalarField::
heatFluxFromCatalyticSurfacesFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict)
{}


Foam::heatFluxFromCatalyticSurfacesFvPatchScalarField::
heatFluxFromCatalyticSurfacesFvPatchScalarField
(
    const heatFluxFromCatalyticSurfacesFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf)
{}


Foam::heatFluxFromCatalyticSurfacesFvPatchScalarField::
heatFluxFromCatalyticSurfacesFvPatchScalarField
(
    const heatFluxFromCatalyticSurfacesFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heatFluxFromCatalyticSurfacesFvPatchScalarField::updateCoeffs()
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

    const psiReactionThermophysicalTransportModel& thermophysicalTransport =
    db().lookupObject<psiReactionThermophysicalTransportModel>
    (
        IOobject::groupName("thermophysicalTransport", IOobject::group(internalField().name()))
    );

    const scalarField Qdotw = surfaceReaction.heatFluxAtCatalyticSurfaces(patchi); // [J/(m2.s)]

    const scalarField rhoDw = thermophysicalTransport.alphaEff(patchi);

  
    const scalarField kappaw = thermo.kappa(patchi); // [J/(m.s.K)]
    volScalarField Cp = thermo.Cp();
    scalarField Cpw = Cp.boundaryField()[patchi];    

    scalarField Qrad(kappaw.size(), 0.0);
    scalar sigma = 5.67e-8; //J/(sec-m2-K4)
    scalar epsilon = 1.0; // assuming a perfect black body

    fvPatchScalarField& Tp =
        const_cast<fvPatchScalarField&>(thermo.T().boundaryField()[patchi]);
    
   forAll(Tp, facei)
   {
       label faceCelli = patch().faceCells()[facei]; 
       Qrad[facei] = sigma*epsilon*
              ( pow(thermo.T().primitiveField()[faceCelli], 4.0) - pow(Tp[facei], 4.0));
   }

    gradient() = (Qdotw-Qrad)/(kappaw);

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::heatFluxFromCatalyticSurfacesFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        heatFluxFromCatalyticSurfacesFvPatchScalarField
    );
}

// ************************************************************************* //
