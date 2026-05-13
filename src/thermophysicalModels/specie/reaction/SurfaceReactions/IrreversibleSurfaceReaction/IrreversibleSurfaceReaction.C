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

#include "IrreversibleSurfaceReaction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class GasThermoType, class SolidThermoType, class ReactionRate>
Foam::IrreversibleSurfaceReaction<GasThermoType, SolidThermoType, ReactionRate>::
IrreversibleSurfaceReaction
(
    const SurfaceReaction<GasThermoType, SolidThermoType>& reaction,
    const ReactionRate& k
)
:
    SurfaceReaction<GasThermoType, SolidThermoType>(reaction),
    k_(k)
{}


template<class GasThermoType, class SolidThermoType, class ReactionRate>
Foam::IrreversibleSurfaceReaction<GasThermoType, SolidThermoType, ReactionRate>::
IrreversibleSurfaceReaction
(
    const speciesTable& gasSpecies,
    const speciesTable& species,    
    const HashPtrTable<GasThermoType>& gasThermoDatabase,
    const HashPtrTable<SolidThermoType>& thermoDatabase,    
    const dictionary& dict
)
:
    SurfaceReaction<GasThermoType, SolidThermoType>
    (
        gasSpecies, 
        species,     
        gasThermoDatabase,
        thermoDatabase,        
        dict
    ),
    k_(species, dict) //check the constructor of ReactionRate later
{}


template<class GasThermoType, class SolidThermoType, class ReactionRate>
Foam::IrreversibleSurfaceReaction<GasThermoType, SolidThermoType, ReactionRate>::
IrreversibleSurfaceReaction
(
    const speciesTable& gasSpecies,
    const speciesTable& species,    
    const HashPtrTable<GasThermoType>& gasThermoDatabase,
    const HashPtrTable<SolidThermoType>& thermoDatabase,    
    const objectRegistry& ob,
    const dictionary& dict
)
:
    SurfaceReaction<GasThermoType, SolidThermoType>
    (
        gasSpecies, 
        species,     
        gasThermoDatabase,
        thermoDatabase,        
        dict
    ),
    k_(species, ob, dict) //check the constructor of ReactionRate later
{}


template<class GasThermoType, class SolidThermoType, class ReactionRate>
Foam::IrreversibleSurfaceReaction<GasThermoType, SolidThermoType, ReactionRate>::
IrreversibleSurfaceReaction
(
    const IrreversibleSurfaceReaction<GasThermoType, SolidThermoType,ReactionRate>& irr,
    const speciesTable& gasSpecies,    
    const speciesTable& species
)
:
    SurfaceReaction<GasThermoType, SolidThermoType>(irr, gasSpecies, species),
    k_(irr.k_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GasThermoType, class SolidThermoType, class ReactionRate>
Foam::scalar Foam::IrreversibleSurfaceReaction<GasThermoType, SolidThermoType, ReactionRate>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    List<scalar> siteDensityList,
    List<label> siteNumberList    
) const
{
    const scalar pi = 3.14159265359;
    const scalar RRjoule = 8.31451; // J/mol-K
    // const scalar RRcal = 1.987316;  // cal/mol-K
        
    scalar stick = 1;
    scalar covModification = 1;
    scalar Tcoeff = 1;
    
    if ( k_.isStick() )
    {
        //- update stick
        scalar coeff1 = 1;
        scalar coeff2 = 1;
        label  k = 0;
        scalar m = 0;
        scalar gamma = 0;
 
        forAll(this->species(), i)
        {
            gamma += c[this->nGasSpecies_+i];
        }

        forAll(this->lhs(), i)
        {
            label j = this->lhs()[i].index;
            scalar order = this->lhs()[i].exponent;
            coeff1 *= pow(siteNumberList[j], order);
            m += this->lhs()[i].stoichCoeff;
        }

        coeff1 /= pow(gamma, m);

        forAll(this->lhsG(), i)
        {
            if (this->lhsG()[i].isgas)
            {
                k = this->lhsG()[i].index;
                coeff2 = sqrt(1000*RRjoule*T/(2*pi*this->gasWk_[k]));
            }

            if (i > 0)
            {
               FatalErrorInFunction 
                   << "The number of gas species in the sticking coefficient "
                   << "reaction cannot be larger than 1." << endl 
                   << "Plase check again your sticking coeffficient reaction" << endl
                   << exit(FatalError);
            }
        }

        stick = coeff1*coeff2;   
    }

    if (k_.isCoverageDependence())
    {
        //- update for covModification
        forAll(k_.coverageModificationCoeffs(), i)
        {  
            scalar etai = k_.coverageModificationCoeffs()[i].second()[0];
            scalar mui  = k_.coverageModificationCoeffs()[i].second()[1];
            scalar epsi = k_.coverageModificationCoeffs()[i].second()[2];
            scalar Zi   = c[this->nGasSpecies_+i]/siteDensityList[i];
            
            //covModification *= pow(10.0, etai*Zi)*pow(Zi, mui)*exp(-epsi*Zi/(RRjoule*T));
            covModification *= pow(10.0, etai*Zi)*pow(Zi, mui)*exp(-epsi*Zi/(T));
        }
    }

    if (k_.isTakahashi())
    {
        //- update coefficient for Takahashi model
        scalar sumCg_ = 0;

        if (this->lhsG().size() < 1 )
        {
            sumCg_ = 1;
        }
        else
        {
            for( label i = 0; i < this->nGasSpecies_;  i++ )
            {
                sumCg_ += c[i];
            }
        }

        Tcoeff = 1.0/sumCg_;
    }

    return stick*covModification*Tcoeff*k_(p, T, c, li);

}


template<class GasThermoType, class SolidThermoType, class ReactionRate>
Foam::scalar Foam::IrreversibleSurfaceReaction<GasThermoType, SolidThermoType, ReactionRate>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return 0;
}


template<class GasThermoType, class SolidThermoType, class ReactionRate>
Foam::scalar Foam::IrreversibleSurfaceReaction<GasThermoType, SolidThermoType, ReactionRate>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return 0;
}


template<class GasThermoType, class SolidThermoType, class ReactionRate>
Foam::scalar Foam::IrreversibleSurfaceReaction<GasThermoType, SolidThermoType, ReactionRate>::dkfdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return k_.ddT(p, T, c, li);
}


template<class GasThermoType, class SolidThermoType, class ReactionRate>
Foam::scalar Foam::IrreversibleSurfaceReaction<GasThermoType, SolidThermoType, ReactionRate>::dkrdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const scalar dkfdT,
    const scalar kr
) const
{
    return 0;
}


template<class GasThermoType, class SolidThermoType, class ReactionRate>
const Foam::List<Foam::Tuple2<Foam::label, Foam::scalarField>>&
Foam::IrreversibleSurfaceReaction<GasThermoType, SolidThermoType, ReactionRate>::coverageModificationCoeffs() const
{
    return k_.coverageModificationCoeffs();
}


template<class GasThermoType, class SolidThermoType, class ReactionRate>
void Foam::IrreversibleSurfaceReaction<GasThermoType, SolidThermoType, ReactionRate>::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcidc
) const
{
    k_.dcidc(p, T, c, li, dcidc);
}


template<class GasThermoType, class SolidThermoType, class ReactionRate>
Foam::scalar Foam::IrreversibleSurfaceReaction<GasThermoType, SolidThermoType, ReactionRate>::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return k_.dcidT(p, T, c, li);
}


template<class GasThermoType, class SolidThermoType, class ReactionRate>
void Foam::IrreversibleSurfaceReaction<GasThermoType, SolidThermoType, ReactionRate>::write
(
    Ostream& os
) const
{
    SurfaceReaction<GasThermoType, SolidThermoType>::write(os);
    k_.write(os);
}


// ************************************************************************* //
