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

#include "DTMMultiComponentMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
void Foam::DTMMultiComponentMixture<ThermoType>::correctMassFractions()
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
Foam::DTMMultiComponentMixture<ThermoType>::DTMMultiComponentMixture
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
    linearityMk_(this->species_.size()),
    epsilonOverKbMk_(this->species_.size()),
    sigmaMk_(this->species_.size()),
    miuiMk_(this->species_.size()),
    polarMk_(this->species_.size()),
    ZrotMk_(this->species_.size()),
    CpCoeffTableMk_(this->species_.size()),

    EPSILONijOVERKB_(this->species_.size()),
    DELTAij_(this->species_.size()),
    Mij_(this->species_.size()),
    SIGMAij_(this->species_.size()),
    mesh_(mesh)
    //
{
    correctMassFractions();

    // precalculation for kinetic theory model
    //for Kinetic model
    forAll(ListW_, i)
    { 
        ListW_[i]           = this->specieThermos()[i].W();
        linearityMk_[i]     = this->specieThermos()[i].linearity();
        epsilonOverKbMk_[i] = this->specieThermos()[i].epsilonOverKb();
        sigmaMk_[i]         = this->specieThermos()[i].sigma();
        miuiMk_[i]          = this->specieThermos()[i].miui();
        polarMk_[i]         = this->specieThermos()[i].polar();
        ZrotMk_[i]          = this->specieThermos()[i].Zrot();        
        CpCoeffTableMk_[i]  = this->specieThermos()[i].CpCoeffTable();
    }

 
    List<scalar> nEPSILONijOVERKB(this->species_.size());
    List<scalar> nDELTAij(this->species_.size());
    List<scalar> nMij(this->species_.size());
    List<scalar> nSIGMAij(this->species_.size());    

    // Dip = dipole moment 
    const scalar DipMin = 1e-20;

    forAll(EPSILONijOVERKB_, i)
    {
        forAll(nEPSILONijOVERKB, j)
        {
            nMij[j] = 
               1/(1/this->specieThermos()[i].W() + 1/this->specieThermos()[j].W());

            nEPSILONijOVERKB[j] = 
               sqrt(this->specieThermos()[i].epsilonOverKb()*this->specieThermos()[j].epsilonOverKb());

            nSIGMAij[j] = 
               0.5*(this->specieThermos()[i].sigma() + this->specieThermos()[j].sigma());

            // calculate coeficient xi
            scalar xi = 1.0; 
            if ((this->specieThermos()[i].miui() < DipMin) && (this->specieThermos()[j].miui() > DipMin)) 
            {
                // miui_j > DipMin > miui_i --> j is polar, i is nonpolar              
                nDELTAij[j] = 0;
                xi = 1.0 + 
                     this->specieThermos()[i].polar()*pow(this->specieThermos()[j].miui(), 2)* 
                     sqrt(this->specieThermos()[j].epsilonOverKb()/this->specieThermos()[i].epsilonOverKb())/
                     (4*pow(this->specieThermos()[i].sigma(), 3)*
                            (
                             1e+19*this->specieThermos()[j].Kb()*this->specieThermos()[j].epsilonOverKb()*
                             pow(this->specieThermos()[j].sigma(), 3)
                            )
                     );
            }
            else if ((this->specieThermos()[i].miui() > DipMin) && (this->specieThermos()[j].miui() < DipMin))
            {
                // miui_j < DipMin < miui_i --> i is polar, j is nonpolar
                nDELTAij[j] = 0;
                xi = 1.0 +
                     this->specieThermos()[j].polar()*pow(this->specieThermos()[i].miui(), 2)*
                     sqrt(this->specieThermos()[i].epsilonOverKb()/this->specieThermos()[j].epsilonOverKb())/
                     (4*pow(this->specieThermos()[j].sigma(), 3)*   
                            (
                             1e+19*this->specieThermos()[i].Kb()*this->specieThermos()[i].epsilonOverKb()*
                             pow(this->specieThermos()[i].sigma(), 3)
                            )
                     ); 
            }
            else 
            {            
                xi = 1.0;
                nDELTAij[j] =
                   0.5*(this->specieThermos()[i].miui()*this->specieThermos()[j].miui())/
                   (nEPSILONijOVERKB[j]*1e+19*this->specieThermos()[j].Kb()*pow(nSIGMAij[j], 3));
            }
           
            nEPSILONijOVERKB[j] = nEPSILONijOVERKB[j]*pow(xi, 2);
            nSIGMAij[j] = nSIGMAij[j]*pow(xi, -1/6);
        }

        EPSILONijOVERKB_[i] = nEPSILONijOVERKB;
        Mij_[i]             = nMij;
        SIGMAij_[i]         = nSIGMAij;
        DELTAij_[i]         = nDELTAij;
    }
    // End of pre-calculation for kinetic model

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::DTMMultiComponentMixture<ThermoType>::cellMixture
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
        Y, X, EPSILONijOVERKB_, DELTAij_, Mij_, SIGMAij_,
        linearityMk_, epsilonOverKbMk_, sigmaMk_, miuiMk_, polarMk_, ZrotMk_, ListW_,
        CpCoeffTableMk_
    );
        
    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::DTMMultiComponentMixture<ThermoType>::patchFaceMixture
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
        Y, X, EPSILONijOVERKB_, DELTAij_, Mij_, SIGMAij_,
        linearityMk_, epsilonOverKbMk_, sigmaMk_, miuiMk_, polarMk_, ZrotMk_, ListW_,
        CpCoeffTableMk_
    );
   
    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::DTMMultiComponentMixture<ThermoType>::cellVolMixture
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
const ThermoType& Foam::DTMMultiComponentMixture<ThermoType>::
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
