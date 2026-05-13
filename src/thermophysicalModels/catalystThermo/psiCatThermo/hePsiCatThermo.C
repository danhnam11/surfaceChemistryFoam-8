/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "hePsiCatThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiCatThermo<BasicPsiThermo, MixtureType>::calculate()
{
    //Info << "using calculate " << endl;
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );
        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];       
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);
                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
  
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::hePsiCatThermo<BasicPsiThermo, MixtureType>::hePsiCatThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heCatThermo<BasicPsiThermo, MixtureType>(mesh, phaseName),
    numberOfSpecies_(MixtureType::numberOfSpecies()), //Nam
    Dimix_(MixtureType::numberOfSpecies()), //
    DimixT_(MixtureType::numberOfSpecies()), //
    Dij_(MixtureType::numberOfSpecies()), //
    specieThermos_(this->specieThermos()),
    mesh_(mesh)
{
    //
    forAll(Dimix_, i)
    {
        Dimix_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    this->phasePropertyName("thermo:Dimix"),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionSet(0, 2, -1, 0, 0)
            )
        );
    }

    forAll(DimixT_, i)
    {
        DimixT_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    this->phasePropertyName("thermo:DimixT"),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionSet(1, -1, -1, 0, 0)
            )
        );
    }
    forAll(Dij_, i)
    {
        Dij_[i] = PtrList<volScalarField>(MixtureType::numberOfSpecies());

        forAll(Dij_[i], j)
        {
            Dij_[i].set
            (
                j,
                new volScalarField
                (
                    IOobject
                    (
                        this->phasePropertyName("thermo:Dij"),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionSet(0, 2, -1, 0, 0)
                )
            );
        }
    }

    calculate();

    // Switch on saving old time
    this->psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::hePsiCatThermo<BasicPsiThermo, MixtureType>::~hePsiCatThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiCatThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


// ************************************************************************* //
