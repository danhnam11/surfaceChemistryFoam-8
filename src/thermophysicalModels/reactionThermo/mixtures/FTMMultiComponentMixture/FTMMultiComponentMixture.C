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

#include "FTMMultiComponentMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
void Foam::FTMMultiComponentMixture<ThermoType>::correctMassFractions()
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
Foam::FTMMultiComponentMixture<ThermoType>::FTMMultiComponentMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicMultiComponentMixture<ThermoType>
    (
        thermoDict,
        mesh,
        phaseName
    ),
    mixture_("mixture", this->specieThermos()[0]),
    mixtureVol_("volMixture", this->specieThermos()[0]),
    // for diffusivity
    ListW_(this->species_.size()),
    muCoeffsMk_(this->species_.size()),
    kappaCoeffsMk_(this->species_.size()),
    DijCoeffsMk_(this->species_.size()),
    mesh_(mesh)
    //
{
    correctMassFractions();

    // precalculation for kinetic theory model
    //for Kinetic model
    forAll(ListW_, i)
    { 
        ListW_[i]           = this->specieThermos()[i].W();
        //muCoeffsMk_[i]      = this->specieThermos()[i].muCoeffs();
        //kappaCoeffsMk_[i]   = this->specieThermos()[i].kappaCoeffs();
        //DijCoeffsMk_[i]     = this->specieThermos()[i].DijCoeffs();
        muCoeffsMk_[i]      = this->specieThermos()[i].coeffs().muCoeffs();
        kappaCoeffsMk_[i]   = this->specieThermos()[i].coeffs().kappaCoeffs();
        DijCoeffsMk_[i]     = this->specieThermos()[i].coeffs().DijCoeffs();
    }

    // End of pre-calculation for kinetic model
    //

    mixture_.updateTRANSFitCoeff
    (
        muCoeffsMk_, kappaCoeffsMk_, DijCoeffsMk_
    );

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::FTMMultiComponentMixture<ThermoType>::cellMixture
(
    const label celli
) const
{
    mixture_ = this->Y_[0][celli]*this->specieThermos()[0];

    for (label n=1; n<this->Y_.size(); n++)
    {
        mixture_ += this->Y_[n][celli]*this->specieThermos()[n];
    }
 
   //- Update coefficients for diffusivity mixture
    //- List of secies mole and mass fraction 
    List<scalar> X(this->Y_.size()); 
    List<scalar> Y(this->Y_.size()); 
    scalar sumXb = 0.0;  
    forAll(X, i)
    {
        sumXb = sumXb + this->Y_[i][celli]/ListW_[i]; 
    }  
    if (sumXb == 0){ sumXb = 1e-30;} 

    forAll(X, i)
    {
        X[i] = (this->Y_[i][celli]/ListW_[i])/sumXb;
        Y[i] = this->Y_[i][celli];
        if(X[i] <= 0) { X[i] = 0; }
        if(Y[i] <= 0) { Y[i] = 0; }
    }

    scalar WmixCorrect = 0.0, sumXcorrected = 0.0;
    forAll(X, i)
    {
        X[i] = X[i] + 1e-40;
        sumXcorrected = sumXcorrected + X[i];
    }
    
    forAll(X, i)
    {
        X[i] = X[i]/sumXcorrected;
        WmixCorrect = WmixCorrect + X[i]*ListW_[i];
    }

    forAll(Y, i)
    {
        Y[i] = X[i]*ListW_[i]/WmixCorrect;
    }
    // Update coefficients for mixture of Kinetic model
    mixture_.updateTRANS
    (
        Y, X, ListW_
    );

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::FTMMultiComponentMixture<ThermoType>::patchFaceMixture
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

   //- Update coefficients for diffusivity mixture
    //- List of secies mole and mass fraction 
    List<scalar> X(this->Y_.size());
    List<scalar> Y(this->Y_.size());
    scalar sumXb = 0.0;
    forAll(X, i)
    {
        sumXb = sumXb + this->Y_[i].boundaryField()[patchi][facei]/ListW_[i];
    }
    if (sumXb == 0){ sumXb = 1e-30;}

    forAll(X, i)
    {
        X[i] = (this->Y_[i].boundaryField()[patchi][facei]/ListW_[i])/sumXb;
        Y[i] = this->Y_[i].boundaryField()[patchi][facei];
        if(X[i] <= 0) { X[i] = 0; } 
        if(Y[i] <= 0) { Y[i] = 0; }
    }

    scalar WmixCorrect = 0.0, sumXcorrected = 0.0;
    forAll(X, i)
    {
        X[i] = X[i] + 1e-40;
        sumXcorrected = sumXcorrected + X[i];
    }

    forAll(X, i)
    {
        X[i] = X[i]/sumXcorrected;
        WmixCorrect = WmixCorrect + X[i]*ListW_[i];
    }

    forAll(Y, i)
    {
        Y[i] = X[i]*ListW_[i]/WmixCorrect;
    }
   
    // Update coefficients for mixture of Kinetic model
    mixture_.updateTRANS
    (
        Y, X, ListW_
    );

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::FTMMultiComponentMixture<ThermoType>::cellVolMixture
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
const ThermoType& Foam::FTMMultiComponentMixture<ThermoType>::
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
