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

#include "hePsiThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiThermo<BasicPsiThermo, MixtureType>::calculate()
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



template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiThermo<BasicPsiThermo, MixtureType>::initialize()
{
    Info << "using Detailed Transport Model without CoTHERM " << endl;
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

    scalarField& CpCells = CpNew_.primitiveFieldRef(); //
    scalarField& kappaCells = kappaNew_.primitiveFieldRef(); //
    scalarField& WmixCells = WmixNew_.primitiveFieldRef(); //

    scalarField& TCellsOld = this->T_.oldTime().primitiveFieldRef();

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

        residualT_.primitiveFieldRef()[celli] = TCells[celli] - TCellsOld[celli]; 

        psiCells[celli]   = mixture_.psi(pCells[celli], TCells[celli]);
        muCells[celli]    = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
        CpCells[celli]    = mixture_.Cp(pCells[celli], TCells[celli]);
        kappaCells[celli] = mixture_.kappa(pCells[celli], TCells[celli]);
        WmixCells[celli]  = mixture_.W();
        
        forAll(Dimix_, i)
        {
            Dimix_[i].primitiveFieldRef()[celli] 
          = mixture_.Dimix(i, pCells[celli], TCells[celli]);

            DimixT_[i].primitiveFieldRef()[celli]
          = mixture_.DimixT(i, pCells[celli], TCells[celli]);

            this->heList_[i].primitiveFieldRef()[celli]
          = specieThermos_[i].HE(pCells[celli], TCells[celli]);
        }
    }

    volScalarField::Boundary& pBf = this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
    volScalarField::Boundary& heBf = this->he().boundaryFieldRef();
    volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
    volScalarField::Boundary& alphaBf = this->alpha_.boundaryFieldRef();
    volScalarField::Boundary& CpBf = CpNew_.boundaryFieldRef();
    volScalarField::Boundary& kappaBf = kappaNew_.boundaryFieldRef();
    volScalarField::Boundary& WmixBf = WmixNew_.boundaryFieldRef();

    // old field 
    volScalarField::Boundary& TBfOld = this->T_.oldTime().boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];       
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        fvPatchScalarField& pWmix = WmixBf[patchi];

        // old field
        fvPatchScalarField& pTOld = TBfOld[patchi]; 

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                residualT_.boundaryFieldRef()[patchi][facei] = pT[facei] - pTOld[facei];

                phe[facei]    = mixture_.HE(pp[facei], pT[facei]);
                ppsi[facei]   = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei]    = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                pCp[facei]    = mixture_.Cp(pp[facei], pT[facei]);
                pkappa[facei] = mixture_.kappa(pp[facei], pT[facei]);
                pWmix[facei]  = mixture_.W();

                forAll(Dimix_, i)
                {
                    Dimix_[i].boundaryFieldRef()[patchi][facei]
                  = mixture_.Dimix(i, pp[facei], pT[facei]);

                    DimixT_[i].boundaryFieldRef()[patchi][facei]
                  = mixture_.DimixT(i, pp[facei], pT[facei]);

                    this->heList_[i].boundaryFieldRef()[patchi][facei]
                  = specieThermos_[i].HE(pp[facei], pT[facei]);
                }

            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);

                residualT_.boundaryFieldRef()[patchi][facei] = pT[facei] - pTOld[facei];

                ppsi[facei]   = mixture_.psi(pp[facei], pT[facei]);
                pmu[facei]    = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                pCp[facei]    = mixture_.Cp(pp[facei], pT[facei]);
                pkappa[facei] = mixture_.kappa(pp[facei], pT[facei]);
                pWmix[facei]  = mixture_.W();

                forAll(Dimix_, i)
                {
                    Dimix_[i].boundaryFieldRef()[patchi][facei]
                  = mixture_.Dimix(i, pp[facei], pT[facei]);

                    DimixT_[i].boundaryFieldRef()[patchi][facei]
                  = mixture_.DimixT(i, pp[facei], pT[facei]);

                    this->heList_[i].boundaryFieldRef()[patchi][facei]
                  = specieThermos_[i].HE(pp[facei], pT[facei]);
                }
                
            }
        }
    }
}



template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiThermo<BasicPsiThermo, MixtureType>::calculateUsingCoTHERMOnlyT()
{
    Info << "using CoTHERM with Only T " << endl;

    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();
    scalarField& CpCells = CpNew_.primitiveFieldRef();
    scalarField& kappaCells = kappaNew_.primitiveFieldRef();
    scalarField& WmixCells = WmixNew_.primitiveFieldRef();

    scalarField& TCellsOld = this->T_.oldTime().primitiveFieldRef();
    
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

        residualT_.primitiveFieldRef()[celli] = TCells[celli] - TCellsOld[celli]; 

        if ( mag(TCells[celli] - TCellsOld[celli] ) < epsilonOnlyT_ )   
        {
            // don't need to calculate thermophysical properties 
            // just copy from old time fields 
            //flagT_.primitiveFieldRef()[celli] = 0.0;
            psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]); 
            muCells[celli] = this->mu_.oldTime()[celli];
            alphaCells[celli] = this->alpha_.oldTime()[celli];
            CpCells[celli] = CpNew_.oldTime()[celli];
            kappaCells[celli] = kappaNew_.oldTime()[celli];
            WmixCells[celli] = mixture_.W();

            forAll(Dimix_, i)
            {
                Dimix_[i].primitiveFieldRef()[celli] 
              = Dimix_[i].oldTime()[celli]; 

                DimixT_[i].primitiveFieldRef()[celli]
              = DimixT_[i].oldTime()[celli];

                this->heList_[i].primitiveFieldRef()[celli]
              = this->heList_[i].oldTime()[celli];
    	    }
        }
        else
        {
            // otherwise, calculate the properties again 
            //flagT_.primitiveFieldRef()[celli] = 1.0;
            psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
            muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
            alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
            CpCells[celli] = mixture_.Cp(pCells[celli], TCells[celli]);
            kappaCells[celli] = mixture_.kappa(pCells[celli], TCells[celli]);
            WmixCells[celli] = mixture_.W();
         
            forAll(Dimix_, i)
            {
                Dimix_[i].primitiveFieldRef()[celli] 
              = mixture_.Dimix(i, pCells[celli], TCells[celli]);

                DimixT_[i].primitiveFieldRef()[celli]
              = mixture_.DimixT(i, pCells[celli], TCells[celli]);

                this->heList_[i].primitiveFieldRef()[celli]
              = specieThermos_[i].HE(pCells[celli], TCells[celli]);
    	    }
        }

    }

    volScalarField::Boundary& pBf = this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
    volScalarField::Boundary& heBf = this->he().boundaryFieldRef();
    volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
    volScalarField::Boundary& alphaBf = this->alpha_.boundaryFieldRef();
    volScalarField::Boundary& CpBf = CpNew_.boundaryFieldRef();
    volScalarField::Boundary& kappaBf = kappaNew_.boundaryFieldRef();
    volScalarField::Boundary& WmixBf = WmixNew_.boundaryFieldRef();
    // old fields
    volScalarField::Boundary& TBfOld = this->T_.oldTime().boundaryFieldRef();
    volScalarField::Boundary& muBfOld = this->mu_.oldTime().boundaryFieldRef();
    volScalarField::Boundary& alphaBfOld = this->alpha_.oldTime().boundaryFieldRef();
    volScalarField::Boundary& CpBfOld = CpNew_.oldTime().boundaryFieldRef();
    volScalarField::Boundary& kappaBfOld = kappaNew_.oldTime().boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];       
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        fvPatchScalarField& pWmix = WmixBf[patchi];

        // old fields
        fvPatchScalarField& pTOld = TBfOld[patchi]; 
        fvPatchScalarField& pmuOld = muBfOld[patchi];
        fvPatchScalarField& palphaOld = alphaBfOld[patchi];
        fvPatchScalarField& pCpOld = CpBfOld[patchi];
        fvPatchScalarField& pkappaOld = kappaBfOld[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);
   
                residualT_.boundaryFieldRef()[patchi][facei] = pT[facei] - pTOld[facei];

                if ( mag(pT[facei] - pTOld[facei] ) < epsilonOnlyT_ )
                {
                    // If the temperature is unchanged,  
                    // don't need to recalculate transport properties.
                    //flagT_.boundaryFieldRef()[patchi][facei] = 0.0;
                    phe[facei] = mixture_.HE(pp[facei], pT[facei]);
                    ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                    pWmix[facei] = mixture_.W();
                    pmu[facei] = pmuOld[facei];
                    palpha[facei] = palphaOld[facei];
                    pCp[facei] = pCpOld[facei];
                    pkappa[facei] = pkappaOld[facei];

                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = Dimix_[i].oldTime().boundaryFieldRef()[patchi][facei];

                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = DimixT_[i].oldTime().boundaryFieldRef()[patchi][facei];

                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = this->heList_[i].oldTime().boundaryFieldRef()[patchi][facei];
                    }
                }
                else 
                {
                    // otherwise, calculate the properties again 
                    //flagT_.boundaryFieldRef()[patchi][facei] = 1.0;
                    phe[facei] = mixture_.HE(pp[facei], pT[facei]);
                    ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                    pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                    palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                    pCp[facei] = mixture_.Cp(pp[facei], pT[facei]);
                    pkappa[facei] = mixture_.kappa(pp[facei], pT[facei]);
                    pWmix[facei] = mixture_.W();

                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = mixture_.Dimix(i, pp[facei], pT[facei]);

                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = mixture_.DimixT(i, pp[facei], pT[facei]);

                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = specieThermos_[i].HE(pp[facei], pT[facei]);
                    }
                }


            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);

                residualT_.boundaryFieldRef()[patchi][facei] = pT[facei] - pTOld[facei];
                
                if ( mag(pT[facei] - pTOld[facei] ) < epsilonOnlyT_ )
                {
                    // If the temperature is unchanged,  
                    // don't need to recalculate transport properties. 
                    //flagT_.boundaryFieldRef()[patchi][facei] = 0.0;
                    ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                    pWmix[facei] = mixture_.W();
                    pmu[facei] = pmuOld[facei];
                    palpha[facei] = palphaOld[facei];
                    pCp[facei] = pCpOld[facei];
                    pkappa[facei] = pkappaOld[facei];

                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = Dimix_[i].oldTime().boundaryFieldRef()[patchi][facei];

                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = DimixT_[i].oldTime().boundaryFieldRef()[patchi][facei];

                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = this->heList_[i].oldTime().boundaryFieldRef()[patchi][facei];

                    }
                }
                else 
                {
                    // otherwise, calculate the properties again 
                    //flagT_.boundaryFieldRef()[patchi][facei] = 1.0;
                    ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                    pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                    palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                    pCp[facei] = mixture_.Cp(pp[facei], pT[facei]);
                    pkappa[facei] = mixture_.kappa(pp[facei], pT[facei]);
                    pWmix[facei] = mixture_.W();

                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = mixture_.Dimix(i, pp[facei], pT[facei]);

                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = mixture_.DimixT(i, pp[facei], pT[facei]);

                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = specieThermos_[i].HE(pp[facei], pT[facei]);
                    }
                }

    
            }
        }
    }
    //Info << "Nam: finish calculating properties " << endl;
}



template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiThermo<BasicPsiThermo, MixtureType>::calculateUsingCoTHERM()
{
    Info << "using CoTHERM " << endl;

    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();
    scalarField& CpCells = CpNew_.primitiveFieldRef();
    scalarField& kappaCells = kappaNew_.primitiveFieldRef();
    scalarField& WmixCells = WmixNew_.primitiveFieldRef();

    // For CoTHERMStepCount
    scalarField& CountCells = coTHERMStepCount_.primitiveFieldRef();

    scalarField& TCellsOld = this->T_.oldTime().primitiveFieldRef();
    scalarField& PCellsOld = this->p_.oldTime().primitiveFieldRef();

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

        residualT_.primitiveFieldRef()[celli] = TCells[celli] - TCellsOld[celli]; 
        scalar coDeltaT = mag(TCells[celli] - TCellsOld[celli]);
        scalar coDeltaP = mag(pCells[celli] - PCellsOld[celli]);

        const bool exceedCount = (CountCells[celli] >= maxCoTHERMStepCount_);

        if 
        ( 
            (coDeltaT <= epsilonT_) && 
            (this->flagSpecies_[celli] < 1.0 ) &&
            (coDeltaP <= epsilonP_) &&
            (!exceedCount)             
        )
        {
            // start counter
            CountCells[celli] += 1.0;

            // don't need to re-calculate thermophysical properties
            // just copy from old time 
            flagCell_.primitiveFieldRef()[celli] = 0.0;
            psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]); 
            muCells[celli] = this->mu_.oldTime()[celli];
            alphaCells[celli] = this->alpha_.oldTime()[celli];
            CpCells[celli] = CpNew_.oldTime()[celli];
            kappaCells[celli] = kappaNew_.oldTime()[celli];
            WmixCells[celli] = mixture_.W();
            forAll(Dimix_, i)
            {
                Dimix_[i].primitiveFieldRef()[celli] 
              = Dimix_[i].oldTime()[celli]; 

                DimixT_[i].primitiveFieldRef()[celli]
              = DimixT_[i].oldTime()[celli];

                this->heList_[i].primitiveFieldRef()[celli]
              = this->heList_[i].oldTime()[celli];
            }

        }
        else if 
        (
            (coDeltaT <= epsilonT_) && 
            (this->flagSpecies_[celli] < 1.0 ) &&
            !(coDeltaP <= epsilonP_) &&
            (!exceedCount)            
        )
        {
            // reset counter
            CountCells[celli] = 0.0;

            // only re-calculate Dimix
            // for other thermophysical properties, copy from old time
            flagCell_.primitiveFieldRef()[celli] = 0.5;
            psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]); 
            muCells[celli] = this->mu_.oldTime()[celli];
            alphaCells[celli] = this->alpha_.oldTime()[celli];
            CpCells[celli] = CpNew_.oldTime()[celli];
            kappaCells[celli] = kappaNew_.oldTime()[celli];
            WmixCells[celli] = mixture_.W();
            forAll(Dimix_, i)
            {
                Dimix_[i].primitiveFieldRef()[celli] 
              = mixture_.Dimix(i, pCells[celli], TCells[celli]);

                DimixT_[i].primitiveFieldRef()[celli]
              = DimixT_[i].oldTime()[celli];

                this->heList_[i].primitiveFieldRef()[celli]
              = this->heList_[i].oldTime()[celli];
            }

        }
        else
        {
            // reset counter
            CountCells[celli] = 0.0;

            // otherwise, re-calculate all thermophysical properties 
            flagCell_.primitiveFieldRef()[celli] = 1.0;
            psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
            muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
            alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
            CpCells[celli] = mixture_.Cp(pCells[celli], TCells[celli]);
            kappaCells[celli] = mixture_.kappa(pCells[celli], TCells[celli]);
            WmixCells[celli] = mixture_.W();

            forAll(Dimix_, i)
            {
                Dimix_[i].primitiveFieldRef()[celli] 
              = mixture_.Dimix(i, pCells[celli], TCells[celli]);

                DimixT_[i].primitiveFieldRef()[celli]
              = mixture_.DimixT(i, pCells[celli], TCells[celli]);

                this->heList_[i].primitiveFieldRef()[celli]
              = specieThermos_[i].HE(pCells[celli], TCells[celli]);
            }

        }

    }

    volScalarField::Boundary& pBf = this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
    volScalarField::Boundary& heBf = this->he().boundaryFieldRef();
    volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
    volScalarField::Boundary& alphaBf = this->alpha_.boundaryFieldRef();
    volScalarField::Boundary& CpBf = CpNew_.boundaryFieldRef();
    volScalarField::Boundary& kappaBf = kappaNew_.boundaryFieldRef();
    volScalarField::Boundary& WmixBf = WmixNew_.boundaryFieldRef();

    volScalarField::Boundary& CountBf  = coTHERMStepCount_.boundaryFieldRef();

    // old fields
    volScalarField::Boundary& TBfOld =
        this->T_.oldTime().boundaryFieldRef();

    volScalarField::Boundary& pBfOld =
        this->p_.oldTime().boundaryFieldRef();

    volScalarField::Boundary& muBfOld =
        this->mu_.oldTime().boundaryFieldRef();

    volScalarField::Boundary& alphaBfOld =
        this->alpha_.oldTime().boundaryFieldRef();

    volScalarField::Boundary& CpBfOld =
        CpNew_.oldTime().boundaryFieldRef();

    volScalarField::Boundary& kappaBfOld =
        kappaNew_.oldTime().boundaryFieldRef();


    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];       
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        fvPatchScalarField& pWmix = WmixBf[patchi];
        fvPatchScalarField& pCount  = CountBf[patchi];

        // old fields
        fvPatchScalarField& pTOld = TBfOld[patchi]; 
        fvPatchScalarField& ppOld = pBfOld[patchi];         
        fvPatchScalarField& pmuOld = muBfOld[patchi];
        fvPatchScalarField& palphaOld = alphaBfOld[patchi];
        fvPatchScalarField& pCpOld = CpBfOld[patchi];
        fvPatchScalarField& pkappaOld = kappaBfOld[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                residualT_.boundaryFieldRef()[patchi][facei] = pT[facei] - pTOld[facei];
                scalar coDeltaTp = mag(pT[facei] - pTOld[facei]);
                scalar coDeltaPp = mag(pp[facei] - ppOld[facei]);

                const bool exceedCountP = (pCount[facei] >= maxCoTHERMStepCount_);

                if 
                ( 
                    (coDeltaTp <= epsilonT_) && 
                    ( this->flagSpecies_.boundaryFieldRef()[patchi][facei] < 1.0) && 
                    (coDeltaPp <= epsilonP_) &&
                    (!exceedCountP)                                       
                )
                {
                    // start counter
                    pCount[facei] += 1.0;

                    // don't need to re-calculate thermophysical properties
                    // just copy from old time
                    flagCell_.boundaryFieldRef()[patchi][facei] = 0.0;
                    phe[facei] = mixture_.HE(pp[facei], pT[facei]);
                    ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                    pWmix[facei] = mixture_.W();
                    pmu[facei] = pmuOld[facei];
                    palpha[facei] = palphaOld[facei];
                    pCp[facei] = pCpOld[facei];
                    pkappa[facei] = pkappaOld[facei];

                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = Dimix_[i].oldTime().boundaryFieldRef()[patchi][facei];

                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = DimixT_[i].oldTime().boundaryFieldRef()[patchi][facei];

                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = this->heList_[i].oldTime().boundaryFieldRef()[patchi][facei];
                    }
                }
                else if 
                (
                    (coDeltaTp <= epsilonT_) && 
                    ( this->flagSpecies_.boundaryFieldRef()[patchi][facei] < 1.0) && 
                    !(coDeltaPp <= epsilonP_) &&
                    (!exceedCountP)                     
                )
                {
                    // reset counter
                    pCount[facei] = 0.0;

                    // only re-calculate Dimix
                    // for other thermophysical properties, copy from old time
                    flagCell_.boundaryFieldRef()[patchi][facei] = 0.5;
                    phe[facei] = mixture_.HE(pp[facei], pT[facei]);
                    ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                    pWmix[facei] = mixture_.W();
                    pmu[facei] = pmuOld[facei];
                    palpha[facei] = palphaOld[facei];
                    pCp[facei] = pCpOld[facei];
                    pkappa[facei] = pkappaOld[facei];

                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = mixture_.Dimix(i, pp[facei], pT[facei]);

                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = DimixT_[i].oldTime().boundaryFieldRef()[patchi][facei];

                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = this->heList_[i].oldTime().boundaryFieldRef()[patchi][facei];
                    }
                }
                else
                {
                    // reset counter
                    pCount[facei] = 0.0;

                    // otherwise, re-calculate all thermophysical properties 
                    flagCell_.boundaryFieldRef()[patchi][facei] = 1.0;
                    phe[facei] = mixture_.HE(pp[facei], pT[facei]);
                    ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                    pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                    palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                    pCp[facei] = mixture_.Cp(pp[facei], pT[facei]);
                    pkappa[facei] = mixture_.kappa(pp[facei], pT[facei]);
                    pWmix[facei] = mixture_.W();

                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = mixture_.Dimix(i, pp[facei], pT[facei]);

                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = mixture_.DimixT(i, pp[facei], pT[facei]);

                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = specieThermos_[i].HE(pp[facei], pT[facei]);
                    }
                }


            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);

                residualT_.boundaryFieldRef()[patchi][facei] = pT[facei] - pTOld[facei];
                scalar coDeltaTp = mag(pT[facei] - pTOld[facei]);
                scalar coDeltaPp = mag(pp[facei] - ppOld[facei]);

                const bool exceedCountP = (pCount[facei] >= maxCoTHERMStepCount_);

                if 
                ( 
                    (coDeltaTp <= epsilonT_) && 
                    ( this->flagSpecies_.boundaryFieldRef()[patchi][facei] < 1.0) &&
                    (coDeltaPp <= epsilonP_) &&
                    (!exceedCountP)                    
                )
                {
                    // start counter
                    pCount[facei] += 1.0;

                    // don't need to re-calculate thermophysical properties
                    // just copy from old time 
                    flagCell_.boundaryFieldRef()[patchi][facei] = 0.0;
                    ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                    pWmix[facei] = mixture_.W();
                    pmu[facei] = pmuOld[facei];
                    palpha[facei] = palphaOld[facei];
                    pCp[facei] = pCpOld[facei];
                    pkappa[facei] = pkappaOld[facei];

                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = Dimix_[i].oldTime().boundaryFieldRef()[patchi][facei];

                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = DimixT_[i].oldTime().boundaryFieldRef()[patchi][facei];

                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = this->heList_[i].oldTime().boundaryFieldRef()[patchi][facei];

                    }
                }
                else if 
                (
                    (coDeltaTp <= epsilonT_) && 
                    ( this->flagSpecies_.boundaryFieldRef()[patchi][facei] < 1.0) &&
                    !(coDeltaPp <= epsilonP_) &&
                    (!exceedCountP)                     
                ) 
                {
                    // reset counter
                    pCount[facei] = 0.0;

                    // only re-calculate Dimix
                    // for other thermophysical properties, copy from old time
                    flagCell_.boundaryFieldRef()[patchi][facei] = 0.5;
                    ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                    pWmix[facei] = mixture_.W();
                    pmu[facei] = pmuOld[facei];
                    palpha[facei] = palphaOld[facei];
                    pCp[facei] = pCpOld[facei];
                    pkappa[facei] = pkappaOld[facei];

                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = mixture_.Dimix(i, pp[facei], pT[facei]);

                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = DimixT_[i].oldTime().boundaryFieldRef()[patchi][facei];

                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = this->heList_[i].oldTime().boundaryFieldRef()[patchi][facei];

                    }

                }
                else
                {
                    // reset counter
                    pCount[facei] = 0.0;
                    
                    // otherwise, re-calculate all thermophysical properties 
                    flagCell_.boundaryFieldRef()[patchi][facei] = 1.0;
                    ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                    pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                    palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                    pCp[facei] = mixture_.Cp(pp[facei], pT[facei]);
                    pkappa[facei] = mixture_.kappa(pp[facei], pT[facei]);
                    pWmix[facei] = mixture_.W();

                    forAll(Dimix_, i)
                    {
                        Dimix_[i].boundaryFieldRef()[patchi][facei]
                      = mixture_.Dimix(i, pp[facei], pT[facei]);

                        DimixT_[i].boundaryFieldRef()[patchi][facei]
                      = mixture_.DimixT(i, pp[facei], pT[facei]);

                        this->heList_[i].boundaryFieldRef()[patchi][facei]
                      = specieThermos_[i].HE(pp[facei], pT[facei]);
                    }
                }

            }
        }
    }
    //Info << "Nam: finish calculating properties " << endl;
}

template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiThermo<BasicPsiThermo, MixtureType>::calculateTransportPreProcessing()
{
    const scalarField& hCells = this->he_;
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();
    scalarField& kappaCells = kappaNew_.primitiveFieldRef();

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
        kappaCells[celli] = mixture_.kappa(pCells[celli], TCells[celli]);

        forAll(Dimix_, i)
        {
            Dimix_[i].primitiveFieldRef()[celli]
          = mixture_.Dimix(i, pCells[celli], TCells[celli]);
        }
        forAll(Dij_, i)
        {
            forAll(Dij_[i],j)
            {
                Dij_[i][j].primitiveFieldRef()[celli]
              = mixture_.Dij(i, j, pCells[celli], TCells[celli]);
            }
        }

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

    volScalarField::Boundary& kappaBf = kappaNew_.boundaryFieldRef();


    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];

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
                pkappa[facei] = mixture_.kappa(pp[facei], pT[facei]);

                forAll(Dimix_, i)
                {
                    Dimix_[i].boundaryFieldRef()[patchi][facei]
                  = mixture_.Dimix(i, pp[facei], pT[facei]);
                }
                forAll(Dij_, i)
                {
                    forAll(Dij_[i], j)
                    {
                        Dij_[i][j].boundaryFieldRef()[patchi][facei]
                      = mixture_.Dij(i, j, pp[facei], pT[facei]);
                    }
                }

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
                pkappa[facei] = mixture_.kappa(pp[facei], pT[facei]);


                forAll(Dimix_, i)
                {
                    Dimix_[i].boundaryFieldRef()[patchi][facei]
                  = mixture_.Dimix(i, pp[facei], pT[facei]);
                }

                forAll(Dij_, i)
                {
                    forAll(Dij_[i], j)
                    {
                        Dij_[i][j].boundaryFieldRef()[patchi][facei]
                      = mixture_.Dij(i, j, pp[facei], pT[facei]);

                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::hePsiThermo<BasicPsiThermo, MixtureType>::hePsiThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName),
    numberOfSpecies_(MixtureType::numberOfSpecies()), //Nam
    Dimix_(MixtureType::numberOfSpecies()), //
    DimixT_(MixtureType::numberOfSpecies()), //
    Dij_(MixtureType::numberOfSpecies()), //
    CpNew_
    (
        IOobject
        (
            "thermo:CpNew",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("CpNew", dimensionSet(0, 2, -2, -1, 0), 0.)
    ),

    kappaNew_
    (
        IOobject
        (
            "kappa_new",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa_new", dimensionSet(1, 1, -3, -1, 0), 0.)
    ),

    WmixNew_
    (
        IOobject
        (
            "Wmix",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Wmix", dimensionSet(1, 0, 0, 0, -1), 0.)
    ),
    specieThermos_(this->specieThermos()),
    mesh_(mesh),

    thermoDict_
    (
        IOdictionary
        (
            IOobject
            (
                this->dictName,
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    ),
    usingDetailedTransportModel_
    (
        thermoDict_.lookupOrDefault("usingDetailedTransportModel", false)
    ),
    usingPreProcessingFTM_
    (
	thermoDict_.lookupOrDefault("usingPreProcessingFTM", false)
    ),
    usingCoTHERM_
    (
        thermoDict_.lookupOrDefault("usingCoTHERM", false)
    ),
    usingCoTHERMOnlyT_
    (
        thermoDict_.lookupOrDefault("usingCoTHERMOnlyT", false)
    ),    
     epsilonOnlyT_
    (
        thermoDict_.lookupOrDefault("epsilonOnlyT", 0.01)
    ),   
    epsilonT_
    (
        thermoDict_.lookupOrDefault("epsilonT", 0.1)
    ),
    epsilonP_
    (
        thermoDict_.lookupOrDefault("epsilonP", 100)
    ),    
    residualT_
    (
        IOobject
        (
            "residual_T",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("residual_T", dimensionSet(0, 0, 0, 1, 0), 0.)
    ),
    residualP_
    (
        IOobject
        (
            "residual_P",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("residual_P", dimensionSet(1, -1, -2, 0, 0), 0.)
    ),    
    flagCell_
    (
        IOobject
        (
            "flagCell",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("flagCell", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    coTHERMStepCount_
    (
        IOobject
        (
            "coTHERMStepCount",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("coTHERMStepCount", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),
    maxCoTHERMStepCount_
    (
        thermoDict_.lookupOrDefault("maxCoTHERMStepCount", 100)
    ),
    mode_(0)    
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

    //Info << "usingDetailedTransportModel = " << this->usingDetailedTransportModel_ << endl;
    //Info << "usingPreProcessingFTM = " << this->usingPreProcessingFTM_ << endl;
    //Info << "usingCoTHERM_ = " << this->usingCoTHERM_ << endl;

    // check and setup coTHERM mode - Nam
    if 
    (
        usingDetailedTransportModel_ && 
        !usingCoTHERM_ && 
        !usingCoTHERMOnlyT_ &&
        !usingPreProcessingFTM_       
    )
    {
        //modeName_ = "DetailedModels";        
        mode_ = 1;
    }
    else if 
    (
        usingDetailedTransportModel_ && 
        usingCoTHERM_ && 
        !usingCoTHERMOnlyT_ &&
        !usingPreProcessingFTM_       
    )
    {
        //modeName_ = "coTHERM";        
        mode_ = 2;
    }
    else if
    (
        usingDetailedTransportModel_ && 
        !usingCoTHERM_ && 
        !usingCoTHERMOnlyT_ &&
        usingPreProcessingFTM_        
    )
    {
        //modeName_ = "preprocessingCoTHERM";
        mode_ = 3;
    }
    else if 
    (
        usingDetailedTransportModel_ && 
        !usingCoTHERM_ && 
        usingCoTHERMOnlyT_ && 
        !usingPreProcessingFTM_       
    )
    {
        //modeName_ = "CoTHERMonlyT";
        mode_ = 4;
    }
    else
    {
        //modeName_ = "originalOpenFOAM";
        mode_ = 0;
    }
    // Nam 

    
    initialize();

    // Switch on saving old time
    this->psi_.oldTime();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::hePsiThermo<BasicPsiThermo, MixtureType>::~hePsiThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    // Nam
    switch(mode_)
    {
        case 1 : 
            initialize();            
            break; 
            
        case 2 : 
            calculateUsingCoTHERM();
            break;      
            
        case 3 : 
            calculateTransportPreProcessing();
            break; 

        case 4 : 
            calculateUsingCoTHERMOnlyT();
            break; 
            
        default:
            calculate();
            break; 
    }


    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


template<class BasicPsiThermo, class MixtureType>
void Foam::hePsiThermo<BasicPsiThermo, MixtureType>::updateOldField
(
    volScalarField Told,
    PtrList<volScalarField> Dimix
)
{
    //Info << "Updating old Field such as: T, mu, kappa, Cp, Dimix" << endl;
    this->Told_ = Told;
    //DimixOld_ = Dimix;
}

// ************************************************************************* //
