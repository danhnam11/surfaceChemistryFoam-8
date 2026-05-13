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

#include "StandardSurfaceChemistryModel.H"
#include "basicMultiComponentMixture.H"
#include "basicMultiComponentCatMixture.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "fvmSup.H" //
#include "localEulerDdtScheme.H" //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
StandardSurfaceChemistryModel
(
    GasReactionThermo& thermo,
    SolidReactionThermo& catThermo    
)
:
    BasicSurfaceChemistryModel<GasReactionThermo, SolidReactionThermo>(thermo, catThermo),
    ODESystem(),
    integrateSurfaceReactionRate_(this->lookupOrDefault("integrateSurfaceReactionRate", true)),
    usingMaxIntegrationTime_(this->lookupOrDefault("usingMaxIntegrationTime", false)),
    maxIntegrationTime_(this->lookupOrDefault("maxIntegrationTime", 1.0)),
    Y_(this->thermo().composition().Y()),
    Ys_(this->catThermo().composition().Y()),    
    gasSpecieThermos_
    (
        dynamic_cast<const basicMultiComponentMixture<GasThermoType>&>
            (this->thermo()).specieThermos()         
    ),
    specieThermos_
    (
        dynamic_cast<const basicMultiComponentCatMixture<SolidThermoType>&>
            (this->catThermo()).specieThermos()        
    ),    
    reactions_
    (
        dynamic_cast<const basicMultiComponentMixture<GasThermoType>&>
        (
            this->thermo()     
        ).species(),
        dynamic_cast<const basicMultiComponentCatMixture<SolidThermoType>&>
        (
            this->catThermo()         
        ).species(),
        gasSpecieThermos_,        
        specieThermos_,
        this->mesh(),
        *this
        //this->subDict("surfaceReactions")        
    ),
    nGasSpecie_(Y_.size()),
    nSolidSpecie_(Ys_.size()),
    nSpecie_(nGasSpecie_ + nSolidSpecie_),
    nReaction_(reactions_.size()),
    Treact_
    (
        BasicSurfaceChemistryModel<GasReactionThermo, SolidReactionThermo>::
        template lookupOrDefault<scalar>
        (
            "Treact",
            0
        )
    ),
    RRg_(nGasSpecie_),
    RRs_(nSolidSpecie_),
    sumRRg_
    (
        IOobject
        (
            "sumRRg",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimMass/dimVolume/dimTime, 0)
    ),        
    c_(nSpecie_+2),
    dcdt_(nSpecie_+2),
    catalystToGeometricAreaRatio_(this->lookupOrDefault("catalystToGeometricAreaRatio", 1.0)),    
    effectivenessFactor_(this->lookupOrDefault("effectivenessFactor", 1.0)),            
    FcatGeo_(catalystToGeometricAreaRatio_*effectivenessFactor_),    
    RRgwall_(nGasSpecie_),
    RRswall_(nSolidSpecie_),
    sumRRgwall_
    (
        this->mesh().boundary(),
        sumRRg_,
        calculatedFvPatchScalarField::typeName
    ),    
    onlySurfaceCatalyst_(this->lookup("onlySurfaceCatalyst")),
    catalyticSurfaces_(this->subDict("catalystRegions").lookup("catalyticSurfaces")),
    siteDensityList_(nSolidSpecie_),
    siteNumberList_(nSolidSpecie_)
{

    forAll(siteDensityList_, i)
    {
        siteDensityList_[i] = specieThermos_[i].siteDensity();
        siteNumberList_[i]  = specieThermos_[i].siteNumber();
    }

    if (integrateSurfaceReactionRate_)
    {
        Info<< "Using integrated surface reaction rate" << endl;
    }
    else
    {
        Info<< "Using instantaneous surface reaction rate" << endl;
    }

    if (usingMaxIntegrationTime_)
    {
        Info<< "Using maxIntegrationTime for surface reaction rate" << endl;
    }
    else
    {}

    // Create the fields for the chemistry sources of gas spcies
    forAll(RRg_, fieldi)
    {
        RRg_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RRg." + Y_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            )
        );

        RRgwall_.set
        (   
            fieldi,
            new volScalarField::Boundary
            (   
                this->mesh().boundary(),
                RRg_[fieldi],
                calculatedFvPatchScalarField::typeName
            )
        );
        forAll(RRgwall_[fieldi], facei)
        {
            RRgwall_[fieldi][facei] = 0.0;
        }
    }

    forAll(RRs_, fieldi)
    {
        RRs_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RRs." + Ys_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            )
        );

        RRswall_.set
        (
            fieldi,
            new volScalarField::Boundary
            (
                this->mesh().boundary(),
                RRs_[fieldi],
                calculatedFvPatchScalarField::typeName
            )
        );
        forAll(RRswall_[fieldi], facei)
        {
            RRswall_[fieldi][facei] = 0.0;
        }

    }

    forAll(sumRRgwall_, patchi)
    {
        forAll(sumRRgwall_[patchi], facei)
        {
            sumRRgwall_[patchi][facei] = 0.0;
        }
    }

    
    if(onlySurfaceCatalyst_)
    {
        Info << "Using Particle Resolved Model for catalyst problem " << endl;
    }
    else
    {
        Info << "Using Porous Media Model for catalyst problem " << endl;
    }
    

    Info << "\nUsing StandardSurfaceChemistryModel for catalyst: " << nl
         << "    number of gas species = " << nGasSpecie_ << nl
         << "    number of solid species = " << nSolidSpecie_ << nl
         << "    number of surface reactions = " << nReaction_ << nl
         << "    site density = " << specieThermos_[0].siteDensity() 
         << " [kmol/m2]" <<  "\n" << endl;

    Info << "catalystToGeometricAreaRatio = " << catalystToGeometricAreaRatio_ << endl; 
    Info << "effectivenessFactor = " << effectivenessFactor_ << endl; 

    /*
    Info << "\nNam check: FcatGeo_ = " << FcatGeo_ << endl;
    Info << "Nam: check initial value of surface coverage" << endl;
    forAll(Ys_, i)
    {
    
        forAll(catalyticSurfaces_, patchI)
        {
            label patchi = this->mesh().boundary().findPatchID(catalyticSurfaces_[patchI]);

            Info << "at patch " << patchi << endl; 
            Info << "i = " << i << " , Ysi = " << Ys_[i].boundaryField()[patchi] << endl;
        }
    }
    */
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
~StandardSurfaceChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
void Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
omega
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcdt
) const
{

    dcdt = Zero;

    forAll(reactions_, i)
    {
        const SurfaceReaction<GasThermoType, SolidThermoType>& R = reactions_[i];

        R.omega(p, T, c, li, dcdt, siteDensityList_, siteNumberList_);
    }
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::scalar Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
omegaI
(
    const label index,
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    const SurfaceReaction<GasThermoType, SolidThermoType>& R = reactions_[index];
    scalar w = R.omega(p, T, c, li, pf, cf, lRef, pr, cr, rRef, siteDensityList_, siteNumberList_);
    return(w);
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
void Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
derivatives
(
    const scalar time,
    const scalarField& c,
    const label li,
    scalarField& dcdt
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    forAll(c_, i)
    {
        c_[i] = max(c[i], 0);
    }

    omega(p, T, c_, li, dcdt);
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
void Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
jacobian
(
    const scalar t,
    const scalarField& c,
    const label li,
    scalarField& dcdt,
    scalarSquareMatrix& J
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    forAll(c_, i)
    {
        c_[i] = max(c[i], 0);
    }

    J = Zero;
    dcdt = Zero;

    scalar omegaI = 0;
    List<label> dummy;
    forAll(reactions_, ri)
    {
        const SurfaceReaction<GasThermoType, SolidThermoType>& R = reactions_[ri];
        scalar kfwd, kbwd;
        R.dwdc
        (
            p, T, c_, li, J, dcdt, omegaI, kfwd, kbwd, false, dummy, siteDensityList_, siteNumberList_
        );
    }

}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::tmp<Foam::volScalarField>
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
tc() const
{ 
    tmp<volScalarField> ttc
    (
        volScalarField::New
        (
            "tc",
            this->mesh(),
            dimensionedScalar(dimTime, small),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    return ttc;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::tmp<Foam::volScalarField>
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
Qdot() const
{
    tmp<volScalarField> tQdot
    (
        volScalarField::New
        (
            "Qdot",
            this->mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        scalarField& Qdot = tQdot.ref();

        forAll(Y_, i)
        {
            forAll(Qdot, celli)
            {
                const scalar hi = specieThermos_[i].Hf();
                Qdot[celli] -= hi*RRs_[i][celli];
            }
        }
    }

    return tQdot;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
calculateRR
(
    const label ri,
    const label si
) const
{
    tmp<volScalarField::Internal> tRR
    (
        volScalarField::Internal::New
        (
            "RR",
            this->mesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );

    volScalarField::Internal& RR = tRR.ref();

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    scalar pf, cf, pr, cr;
    label lRef, rRef;

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = p[celli];

        for (label i=0; i<nSolidSpecie_; i++)
        {
            const scalar Yi = Ys_[i][celli];
            c_[i] = rhoi*Yi/specieThermos_[i].W();
        }

        const SurfaceReaction<GasThermoType, SolidThermoType>& R = reactions_[ri];
        const scalar omegai = R.omega
        (
            pi, Ti, c_, celli, pf, cf, lRef, pr, cr, rRef, siteDensityList_, siteNumberList_
        );

        forAll(R.lhs(), s)
        {
            if (si == R.lhs()[s].index)
            {
                RR[celli] -= R.lhs()[s].stoichCoeff*omegai;
            }
        }

        forAll(R.rhs(), s)
        {
            if (si == R.rhs()[s].index)
            {
                RR[celli] += R.rhs()[s].stoichCoeff*omegai;
            }
        }

        RR[celli] *= specieThermos_[si].W();
    }

    return tRR;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
void Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = p[celli];

        for (label i=0; i<nGasSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c_[i] = rhoi*Yi/gasSpecieThermos_[i].W();
        }

        for (label i=0; i<nSolidSpecie_; i++)
        {
            const scalar Ysi = Ys_[i][celli];
            c_[nGasSpecie_+i] = 
                Ysi*specieThermos_[i].siteDensity()/specieThermos_[i].siteNumber();
        }

        omega(pi, Ti, c_, celli, dcdt_);

        for (label i=0; i<nGasSpecie_; i++)
        {
            RRg_[i][celli] = dcdt_[i]*gasSpecieThermos_[i].W();
        }

        for (label i=0; i<nSolidSpecie_; i++)
        {
            RRs_[i][celli] = dcdt_[nGasSpecie_+i]*specieThermos_[i].W();
        }

    }
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
template<class DeltaTType>
Foam::scalar Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
solve
(
    const DeltaTType& deltaT
)
{
    // For porous media catalyst model with surface reaction
    BasicSurfaceChemistryModel<GasReactionThermo, SolidReactionThermo>::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    scalarField c0(nSpecie_);

    forAll(rho, celli)
    {
        scalar Ti = T[celli];

        if (Ti > Treact_)
        {
            const scalar rhoi = rho[celli];
            scalar pi = p[celli];

            // c of gas species in kmol/m3
            for (label i=0; i<nGasSpecie_; i++)
            {
                c_[i] = rhoi*Y_[i][celli]/gasSpecieThermos_[i].W();
                c0[i] = c_[i];
            }

            // c of site species in kmol/m2
            for (label i=0; i<nSolidSpecie_; i++)
            {
                c_[nGasSpecie_+i] = 
                    Ys_[i][celli]*specieThermos_[i].siteDensity()/specieThermos_[i].siteNumber();
                c0[nGasSpecie_+i] = c_[nGasSpecie_+i];
            }


            // Initialise time progress
            scalar timeLeft = deltaT[celli];

            // Calculate the chemical source terms
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                this->solve(pi, Ti, c_, celli, dt, this->deltaTChem_[celli]);
                timeLeft -= dt;
            }

            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            this->deltaTChem_[celli] =
                min(this->deltaTChem_[celli], this->deltaTChemMax_);

            scalar sumMassGas = 0;
            // RR for gas species in kg/(m3-sec)
            for (label i=0; i<nGasSpecie_; i++)
            {
                RRg_[i][celli] =
                    //(c_[i] - c0[i])*gasSpecieThermos_[i].W()/deltaT[celli];
                    FcatGeo_*(c_[i] - c0[i])*gasSpecieThermos_[i].W()/deltaT[celli];

                // calculate total mass of gas-phase species
                sumMassGas += c_[i]*gasSpecieThermos_[i].W();
            }

            // initialize temporary mass fraction of gas-phase at boundary
            List<scalar> YiBtemp(Y_.size());
            forAll(YiBtemp, i)
            {
                YiBtemp[i] = 0;
            }

            for (label i=0; i<nGasSpecie_; i++)
            {
                // calculate the temporary mass fraction of gas-phase species
                YiBtemp[i] = c_[i]*gasSpecieThermos_[i].W()/sumMassGas;
            }

            scalar sumRRg = 0;
            for (label i=0; i<nGasSpecie_; i++)
            {
                // total gas-phase species flux
                sumRRg += RRg_[i][celli];
            }

            // add correction term induced by stefan flow for specie flux
            for (label i=0; i<nGasSpecie_; i++)
            {
                 RRg_[i][celli] -= YiBtemp[i]*sumRRg;
            }

           // update back Ys_ field after solving surface reactions
            scalar sumCS = 0.0;
            for (label i=0; i<nSolidSpecie_; i++)
            {
                Ys_[i].primitiveFieldRef()[celli] = 
                    max
                    (
                        c_[nGasSpecie_+i]/specieThermos_[i].siteDensity()*
                        specieThermos_[i].siteNumber(), 
                        0.0
                    );

                sumCS += Ys_[i].primitiveField()[celli];
            }

            // normalize 
            for (label i=0; i<nSolidSpecie_; i++)
            {
                Ys_[i].primitiveFieldRef()[celli] = Ys_[i].primitiveFieldRef()[celli]/sumCS;
            }

            // the first solid (site species) = 1 - sum of (n-1 species)
            //Ys_[0].primitiveFieldRef()[celli] = 1.0 - sumCS;

            // RR for site species in 1/sec
            for (label i=0; i<nSolidSpecie_; i++)
            {
                RRs_[i][celli] =
                    (c_[nGasSpecie_+i] - c0[nGasSpecie_+i])*
                    specieThermos_[i].siteNumber()/specieThermos_[i].siteDensity()/deltaT[celli];
            }

        }
    }

    return deltaTMin;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::scalar Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    if (onlySurfaceCatalyst_)
    {
        //  for particle resolved catalyst model with surface reactions
        return min
        (
            this->solvePatches<UniformField<scalar>>(UniformField<scalar>(deltaT)),
            2*deltaT
        );
    }
    else
    {
        // for porous media catalyst model with surface reactions
        return min
        (
            this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
            2*deltaT
        );
    }

}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::scalar Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
solve
(
    const scalarField& deltaT
)
{
    if (onlySurfaceCatalyst_)
    {
        // for particle resolved catalyst model with surface reactions
        return this->solvePatches<scalarField>(deltaT); 
    }
    else
    {
        // for porous media catalyst model with surface reactions
        return this->solve<scalarField>(deltaT);
    }
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
reactionRate
(
    const label ri
) const
{
    tmp<volScalarField::Internal> tRR
    (
        volScalarField::Internal::New
        (
            "RR",
            this->mesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );

    volScalarField::Internal& RR = tRR.ref();

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    scalar pf, cf, pr, cr;
    label lRef, rRef;

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = p[celli];

        for (label i=0; i<nGasSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c_[i] = rhoi*Yi/gasSpecieThermos_[i].W();
        }

        const SurfaceReaction<GasThermoType, SolidThermoType>& R = reactions_[ri];
        const scalar omegai = R.omega
        (
            pi, Ti, c_, celli, pf, cf, lRef, pr, cr, rRef, siteDensityList_, siteNumberList_
        );

        RR[celli] = omegai;
    }

    return tRR;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::tmp<Foam::fvScalarMatrix> Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
RFromSCHEM
(
    volScalarField& gasYi
) const 
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(gasYi, dimMass/dimTime));
    fvScalarMatrix& Su = tSu.ref();

    const label specieI = this->thermo().composition().species()[gasYi.member()];

    Su += RRg(specieI);

    return tSu;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::tmp<Foam::volScalarField> Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
QdotFromSCHEM() const 
{
    tmp<volScalarField> tQdot
    (
        volScalarField::New
        ( 
            "Qdot.FromSurf",
            this->mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );
    
    if (this->chemistry_)
    {
        scalarField& Qdot = tQdot.ref();
        
        forAll(Y_, i)
        {
            forAll(Qdot, celli)
            {
                const scalar hi_g = gasSpecieThermos_[i].Hf();
                Qdot[celli] -= hi_g*RRg_[i][celli];
            }
        }
              
           
        forAll(Ys_, i)
        {   
            forAll(Qdot, celli)
            {
                const scalar hi_s = specieThermos_[i].Hf();
                Qdot[celli] -= hi_s*RRs_[i][celli];
            }
        }   

    }    
                
    return tQdot;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::scalarField Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
speciesFluxAtCatalyticSurfaces
(
    const word speciei,
    const label patchi
) const
{
    label i = this->thermo().composition().species()[speciei];
    return RRgwall_[i][patchi];
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::scalarField Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
heatFluxAtCatalyticSurfaces
(
    const label patchi
) const
{
    scalarField tQdot(this->mesh().boundary()[patchi].size(), 0.0);

    if (this->chemistry_)
    {
        forAll(Y_, i)
        {
            forAll(tQdot, facei)
            {
                scalar hf = gasSpecieThermos_[i].Hf();
                tQdot[facei] -= hf*RRgwall_[i][patchi][facei];
            }
        }
        
        forAll(Ys_, i)
        {
            forAll(tQdot, facei)
            {
                scalar hf = specieThermos_[i].Hf();
                tQdot[facei] -= hf*RRswall_[i][patchi][facei];
            }
            
        }
        
    }

    return tQdot;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::scalarField Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
momentumFluxAtCatalyticSurfaces
(
    const label patchi
) const
{
    return sumRRgwall_[patchi];
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
template<class DeltaTType>
Foam::scalar Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
solvePatches
(
    const DeltaTType& deltaT
)
{
    // for particle resolved catalyst model with surface reactions
    BasicSurfaceChemistryModel<GasReactionThermo, SolidReactionThermo>::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const volScalarField& rho = trho();

    scalarField c0(nSpecie_, 0.0);

    // from catalyticSurface_, find the ID of patches
    forAll(catalyticSurfaces_, patchI)
    {
        label patchi = this->mesh().boundary().findPatchID(catalyticSurfaces_[patchI]);

        const fvPatchScalarField& pT = this->thermo().T().boundaryField()[patchi];
        const fvPatchScalarField& pp = this->thermo().p().boundaryField()[patchi];
        const fvPatchScalarField& pRho = rho.boundaryField()[patchi];        

        forAll(pT, facei)
        { 
            label celli = pT.patch().faceCells()[facei];

            scalar Ti = pT[facei];
            scalar pi = pp[facei];
            scalar gasRhoi = pRho[facei]; 

          if (Ti > Treact_)
          {
            forAll(Y_, i)
            {
                // get value from existing boundary fields of gas-phase species
                RRgwall_[i][patchi][facei] = 0.0;
                scalar Yi = Y_[i].boundaryField()[patchi][facei];
                // then transfer it into c_ vector
                c_[i] = gasRhoi*Yi/gasSpecieThermos_[i].W();
                // then transfer it into c0 vector
                c0[i] = c_[i];
            }
           
            for (label i=0;  i < nSolidSpecie_; i++)
            {
                // get value from existing boundary fields of site species
                scalar Ysi = Ys_[i].boundaryField()[patchi][facei];
                // then transfer it into c_ vector
                c_[nGasSpecie_+i] = 
                    Ysi*specieThermos_[i].siteDensity()/specieThermos_[i].siteNumber();
                // then transfer it into c0 vector 
                c0[nGasSpecie_+i] = c_[nGasSpecie_+i];
            }
           
            // Initialize Stefan flow induced species flux
            sumRRgwall_[patchi][facei] = 0.0;

            // Initialise time pr gress
            scalar timeLeft = deltaT[celli];

            // solving surface chemistry
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                this->solve(pi, Ti, c_, celli, dt, this->deltaTChem_[celli]);
                timeLeft -= dt;
            }

            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);
            this->deltaTChem_[celli] = min(this->deltaTChem_[celli], this->deltaTChemMax_);

            scalar sumMassGas = 0;
            for (label i=0; i<nGasSpecie_; i++)
            {
                // calculate individual gas-phase species production rate from 
                // surface reactions 
                RRgwall_[i][patchi][facei] =
                    FcatGeo_*(c_[i] - c0[i])*gasSpecieThermos_[i].W()/deltaT[celli];
                    //(c_[i] - c0[i])*gasSpecieThermos_[i].W()/deltaT[celli];

                // calculate total mass of gas-phase species
                sumMassGas += c_[i]*gasSpecieThermos_[i].W();
            }

            // initialize temporary mass fraction of gas-phase at boundary
            List<scalar> YiBtemp(Y_.size());
            forAll(YiBtemp, i)
            {
                YiBtemp[i] = 0;
            }

            for (label i=0; i<nGasSpecie_; i++)
            {
                // calculate the temporary mass fraction of gas-phase species
                YiBtemp[i] = c_[i]*gasSpecieThermos_[i].W()/sumMassGas;
            }

            scalar sumRRgwall = 0;
            for (label i=0; i<nGasSpecie_; i++)
            {
                // total gas-phase species flux
                sumRRgwall += RRgwall_[i][patchi][facei];
            }

            // Stefan flow induced species flux
            sumRRgwall_[patchi][facei] = sumRRgwall;

            // add correction term induced by stefan flow for specie flux
            for (label i=0; i<nGasSpecie_; i++)
            {
                 RRgwall_[i][patchi][facei] -= YiBtemp[i]*sumRRgwall;
            }

            // update back Ys_ field after solving surface reactions
            scalar sumCS = 0.0;
            for (label i=1; i<nSolidSpecie_; i++)
            {
                Ys_[i].boundaryFieldRef()[patchi][facei] = 
                    max
                    (
                        c_[nGasSpecie_+i]/specieThermos_[i].siteDensity()*
                        specieThermos_[i].siteNumber(), 
                        0.0
                    );

                sumCS += Ys_[i].boundaryFieldRef()[patchi][facei];
            }


            
            Ys_[0].boundaryFieldRef()[patchi][facei] = 1.0 - sumCS;

            for (label i=0; i<nSolidSpecie_; i++)
            {
                RRswall_[i][patchi][facei] =
                    FcatGeo_*(c_[nGasSpecie_+i] - c0[nGasSpecie_+i])/deltaT[celli];
                    //(c_[nGasSpecie_+i] - c0[nGasSpecie_+i])/deltaT[celli];                    
            }
         
          }//if Ti>Treact_ 

        } // for all facei

    }//forAll patches

    return deltaTMin;
}



template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
void Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
solveSurfaceReactions()
{
    if (integrateSurfaceReactionRate_)
    {
        if (fv::localEulerDdt::enabled(this->mesh()))
        {
            const scalarField& rDeltaT = fv::localEulerDdt::localRDeltaT(this->mesh());
            if (usingMaxIntegrationTime_)
            {
                this->solve( min(1.0/rDeltaT, maxIntegrationTime_)() );
            }
            else
            {
                this->solve((1.0/rDeltaT)());
            }
        }
        else
        {
            this->solve(this->mesh().time().deltaTValue());
        }
    }
    else
    {
        calculate();
    }
}


// ************************************************************************* //
