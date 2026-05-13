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

#include "SpecieCatMixture.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::SpecieCatMixture<MixtureType>::volScalarFieldProperty
(
    const word& psiName,
    const dimensionSet& psiDim,
    scalar (MixtureType::thermoType::*psiMethod)
    (
        const scalar,
        const scalar
    ) const,
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    const typename MixtureType::thermoType& thermo =
        this->specieThermo(speciei);

    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName(psiName, T.group()),
            T.mesh(),
            psiDim
        )
    );

    volScalarField& psi = tPsi.ref();

    forAll(p, celli)
    {
        psi[celli] = (thermo.*psiMethod)(p[celli], T[celli]);
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        const fvPatchScalarField& pp = p.boundaryField()[patchi];
        const fvPatchScalarField& pT = T.boundaryField()[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];

        forAll(pp, facei)
        {
            ppsi[facei] = (thermo.*psiMethod)(pp[facei], pT[facei]);
        }
    }

    return tPsi;
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::SpecieCatMixture<MixtureType>::fieldProperty
(
    scalar (MixtureType::thermoType::*psiMethod)
    (
        const scalar,
        const scalar
    ) const,
    const label speciei,
    const scalarField& p,
    const scalarField& T
) const
{
    const typename MixtureType::thermoType& thermo =
        this->specieThermo(speciei);

    tmp<scalarField> tPsi(new scalarField(p.size()));

    scalarField& psi = tPsi.ref();

    forAll(p, facei)
    {
        psi[facei] = (thermo.*psiMethod)(p[facei], T[facei]);
    }

    return tPsi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::SpecieCatMixture<MixtureType>::SpecieCatMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    MixtureType(thermoDict, mesh, phaseName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType>
Foam::scalar Foam::SpecieCatMixture<MixtureType>::Wi(const label speciei) const
{
    return this->specieThermo(speciei).W();
}


template<class MixtureType>
Foam::scalar Foam::SpecieCatMixture<MixtureType>::Hf(const label speciei) const
{
    return this->specieThermo(speciei).Hf();
}


template<class MixtureType>
Foam::scalar Foam::SpecieCatMixture<MixtureType>::rho
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).rho(p, T);
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::SpecieCatMixture<MixtureType>::rho
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "rho",
        dimDensity,
        &MixtureType::thermoType::rho,
        speciei,
        p,
        T
    );
}


template<class MixtureType>
Foam::scalar Foam::SpecieCatMixture<MixtureType>::Cp
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).Cp(p, T);
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::SpecieCatMixture<MixtureType>::Cp
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "Cp",
        dimEnergy/dimMass/dimTemperature,
        &MixtureType::thermoType::Cp,
        speciei,
        p,
        T
    );
}


template<class MixtureType>
Foam::scalar Foam::SpecieCatMixture<MixtureType>::HE
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).HE(p, T);
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::SpecieCatMixture<MixtureType>::HE
(
    const label speciei,
    const scalarField& p,
    const scalarField& T
) const
{
    return fieldProperty(&MixtureType::thermoType::HE, speciei, p, T);
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::SpecieCatMixture<MixtureType>::HE
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "HE",
        dimEnergy/dimMass,
        &MixtureType::thermoType::HE,
        speciei,
        p,
        T
    );
}


template<class MixtureType>
Foam::scalar Foam::SpecieCatMixture<MixtureType>::Hs
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).Hs(p, T);
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::SpecieCatMixture<MixtureType>::Hs
(
    const label speciei,
    const scalarField& p,
    const scalarField& T
) const
{
    return fieldProperty(&MixtureType::thermoType::Hs, speciei, p, T);
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::SpecieCatMixture<MixtureType>::Hs
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "Hs",
        dimEnergy/dimMass,
        &MixtureType::thermoType::Hs,
        speciei,
        p,
        T
    );
}


template<class MixtureType>
Foam::scalar Foam::SpecieCatMixture<MixtureType>::Ha
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).Ha(p, T);
}


template<class MixtureType>
Foam::tmp<Foam::scalarField> Foam::SpecieCatMixture<MixtureType>::Ha
(
    const label speciei,
    const scalarField& p,
    const scalarField& T
) const
{
    return fieldProperty(&MixtureType::thermoType::Ha, speciei, p, T);
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::SpecieCatMixture<MixtureType>::Ha
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "Ha",
        dimEnergy/dimMass,
        &MixtureType::thermoType::Ha,
        speciei,
        p,
        T
    );
}


template<class MixtureType>
Foam::scalar Foam::SpecieCatMixture<MixtureType>::mu
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).mu(p, T);
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::SpecieCatMixture<MixtureType>::mu
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "mu",
        dimMass/dimLength/dimTime,
        &MixtureType::thermoType::mu,
        speciei,
        p,
        T
    );
}


template<class MixtureType>
Foam::scalar Foam::SpecieCatMixture<MixtureType>::kappa
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).kappa(p, T);
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::SpecieCatMixture<MixtureType>::kappa
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "kappa",
        dimPower/dimLength/dimTemperature,
        &MixtureType::thermoType::kappa,
        speciei,
        p,
        T
    );
}


template<class MixtureType>
Foam::scalar Foam::SpecieCatMixture<MixtureType>::alphah
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).alphah(p, T);
}


template<class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::SpecieCatMixture<MixtureType>::alphah
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldProperty
    (
        "alphah",
        dimMass/dimLength/dimTime,
        &MixtureType::thermoType::alphah,
        speciei,
        p,
        T
    );
}

// ************************************************************************* //
