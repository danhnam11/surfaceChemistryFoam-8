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

#include "multiComponentCatMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
void Foam::multiComponentCatMixture<ThermoType>::correctMassFractions()
{
    // Multiplication by 1.0 changes Yt patches to "calculated"
    volScalarField Yt("Yt", 1.0*this->Y_[0]);

    for (label n=1; n<this->Y_.size(); n++)
    {
        Yt += this->Y_[n];
    }

    if (mag(max(Yt).value()) < rootVSmall)
    {
        FatalErrorInFunction
            << "Sum of mass fractions is zero for species " << this->species()
            << exit(FatalError);
    }
        
    forAll(this->Y_, n)
    {
        this->Y_[n] /= Yt;
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multiComponentCatMixture<ThermoType>::multiComponentCatMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicMultiComponentCatMixture<ThermoType>
    (
        thermoDict,
        mesh,
        phaseName
    ),
    mixture_("mixture", this->specieThermos()[0]),
    mixtureVol_("volMixture", this->specieThermos()[0])  
{
    correctMassFractions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::multiComponentCatMixture<ThermoType>::cellMixture
(
    const label celli
) const
{
    mixture_ = this->Y_[0][celli]*this->specieThermos()[0];

    for (label n=1; n<this->Y_.size(); n++)
    {
        mixture_ += this->Y_[n][celli]*this->specieThermos()[n];
    }

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentCatMixture<ThermoType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    mixture_ = this->Y_[0].boundaryField()[patchi][facei]*this->specieThermos()[0];

    for (label n=1; n<this->Y_.size(); n++)
    {
        mixture_ += this->Y_[n].boundaryField()[patchi][facei]*this->specieThermos()[n];
    }

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentCatMixture<ThermoType>::cellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    scalar rhoInv = 0.0;
    forAll(this->specieThermos(), i)
    {
        rhoInv += this->Y_[i][celli]/this->specieThermos()[i].rho(p, T);
    }

    mixtureVol_ =
        this->Y_[0][celli]/this->specieThermos()[0].rho(p, T)/rhoInv*this->specieThermos()[0];

    for (label n=1; n<this->Y_.size(); n++)
    {
        mixtureVol_ +=
            this->Y_[n][celli]/this->specieThermos()[n].rho(p, T)/rhoInv*this->specieThermos()[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponentCatMixture<ThermoType>::
patchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalar rhoInv = 0.0;
    forAll(this->specieThermos(), i)
    {
        rhoInv +=
            this->Y_[i].boundaryField()[patchi][facei]/this->specieThermos()[i].rho(p, T);
    }

    mixtureVol_ =
        this->Y_[0].boundaryField()[patchi][facei]/this->specieThermos()[0].rho(p, T)/rhoInv
      * this->specieThermos()[0];

    for (label n=1; n<this->Y_.size(); n++)
    {
        mixtureVol_ +=
            this->Y_[n].boundaryField()[patchi][facei]/this->specieThermos()[n].rho(p,T)
          / rhoInv*this->specieThermos()[n];
    }

    return mixtureVol_;
}

// ************************************************************************* //
