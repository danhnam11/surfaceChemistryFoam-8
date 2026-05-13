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

#include "SurfaceReaction.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

template<class GasThermoType, class SolidThermoType>
Foam::label Foam::SurfaceReaction<GasThermoType, SolidThermoType>::
nUnNamedSurfaceReactions(0);

template<class GasThermoType, class SolidThermoType>
Foam::scalar Foam::SurfaceReaction<GasThermoType, SolidThermoType>::
TlowDefault(0);

template<class GasThermoType, class SolidThermoType>
Foam::scalar Foam::SurfaceReaction<GasThermoType, SolidThermoType>::
ThighDefault(great);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GasThermoType, class SolidThermoType>
Foam::label Foam::SurfaceReaction<GasThermoType, SolidThermoType>::
getNewSurfaceReactionID()
{
    return nUnNamedSurfaceReactions++;
}


template<class GasThermoType, class SolidThermoType>
void Foam::SurfaceReaction<GasThermoType, SolidThermoType>::setThermo
(
    const HashPtrTable<SolidThermoType>& thermoDatabase
)
{
    typename SolidThermoType::thermoType rhsThermo
    (
        rhs_[0].stoichCoeff
       *(*thermoDatabase[species_[rhs_[0].index]]).W()
       *(*thermoDatabase[species_[rhs_[0].index]])
    );

    for (label i=1; i<rhs_.size(); ++i)
    {
        rhsThermo +=
            rhs_[i].stoichCoeff
           *(*thermoDatabase[species_[rhs_[i].index]]).W()
           *(*thermoDatabase[species_[rhs_[i].index]]);
    }

    typename SolidThermoType::thermoType lhsThermo
    (
        lhs_[0].stoichCoeff
       *(*thermoDatabase[species_[lhs_[0].index]]).W()
       *(*thermoDatabase[species_[lhs_[0].index]])
    );

    for (label i=1; i<lhs_.size(); ++i)
    {
        lhsThermo +=
            lhs_[i].stoichCoeff
           *(*thermoDatabase[species_[lhs_[i].index]]).W()
           *(*thermoDatabase[species_[lhs_[i].index]]);
    }

    SolidThermoType::thermoType::operator=(lhsThermo == rhsThermo);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
/*
//- Construct from components
template<class GasThermoType, class SolidThermoType>
Foam::SurfaceReaction<GasThermoType, SolidThermoType>::SurfaceReaction
(
    const speciesTable& gasSpecies,
    const speciesTable& species,    
    const List<catSpecieCoeffs>& lhsG,
    const List<catSpecieCoeffs>& rhsG,
    const List<catSpecieCoeffs>& lhs,
    const List<catSpecieCoeffs>& rhs,
    const HashPtrTable<GasThermoType>& gasThermoDatabase,        
    const HashPtrTable<SolidThermoType>& thermoDatabase
)
:
    SolidThermoType::thermoType(*thermoDatabase[species[0]]),
    name_("un-named-reaction-" + Foam::name(getNewSurfaceReactionID())),

    Tlow_(TlowDefault),
    Thigh_(ThighDefault),
    lhs_(lhs),
    rhs_(rhs),
    lhsG_(lhsG), 
    rhsG_(rhsG),    
    gasSpecies_(gasSpecies),
    species_(species),    
    nG_(species.size()), 
    nGasSpecies_(species.size()),
    nSolidSpecies_(gasSpecies.size()),
    gasWk_(nGasSpecies_)
    //thermoDatabase_(thermoDatabase)
{
    setThermo(thermoDatabase);
    forAll(gasWk_, i)
    {
        const scalar gWi = gasThermoDatabase[gasSpecies[i]]->W();
        gasWk_[i] = gWi;
    }
}
*/


//- Construct as copy given new speciesTable
template<class GasThermoType, class SolidThermoType>
Foam::SurfaceReaction<GasThermoType, SolidThermoType>::SurfaceReaction
(
    const SurfaceReaction<GasThermoType, SolidThermoType>& r,
    const speciesTable& gasSpecies,    
    const speciesTable& species
)
:
    SolidThermoType::thermoType(r),
    name_(r.name() + "Copy"),
    Tlow_(r.Tlow()),
    Thigh_(r.Thigh()),
    lhs_(r.lhs_),
    rhs_(r.rhs_),
    lhsG_(r.lhsG_),
    rhsG_(r.rhsG_),
    gasSpecies_(gasSpecies),
    species_(species),    
    nG_(species.size()), 
    nGasSpecies_(gasSpecies.size()),
    nSolidSpecies_(species.size()),
    gasWk_(nGasSpecies_)
    // thermoDatabase_(r.thermoDatabase())
{
    //update this later
    forAll(gasWk_, i)
    {
        gasWk_[i] = 1.0;  //temporary set
    }    
}


//- Construct from dictionary
template<class GasThermoType, class SolidThermoType>
Foam::SurfaceReaction<GasThermoType, SolidThermoType>::SurfaceReaction
(
    const speciesTable& gasSpecies,
    const speciesTable& species,    
    const HashPtrTable<GasThermoType>& gasThermoDatabase,
    const HashPtrTable<SolidThermoType>& thermoDatabase,    
    const dictionary& dict
)
:
    SolidThermoType::thermoType(*thermoDatabase[species[0]]),
    name_(dict.dictName()),
    Tlow_(dict.lookupOrDefault<scalar>("Tlow", TlowDefault)),
    Thigh_(dict.lookupOrDefault<scalar>("Thigh", ThighDefault)),
    gasSpecies_(gasSpecies),
    species_(species),
    nG_(gasSpecies_.size()),
    nGasSpecies_(gasSpecies.size()),
    nSolidSpecies_(species.size()),
    gasWk_(nGasSpecies_)        
{
    forAll(gasWk_, i)
    {
        const scalar gWi = gasThermoDatabase[gasSpecies[i]]->W();
        gasWk_[i] = gWi;
    }

    List<catSpecieCoeffs> lhstemp_;
    List<catSpecieCoeffs> rhstemp_;
    
    catSpecieCoeffs::setLRhsCat
    (
        IStringStream(dict.lookup("reaction"))(),
        gasSpecies_,
        species_,
        lhstemp_,
        rhstemp_
    );

    forAll(lhstemp_,i)
    {
        if(lhstemp_[i].isgas == true) 
        {
            lhsG_.append(lhstemp_[i]);
        }
        else
        {
            lhs_.append(lhstemp_[i]);
        }
    }
    forAll(rhstemp_,i)
    {
        if(rhstemp_[i].isgas == true) 
        {
            rhsG_.append(rhstemp_[i]);
        }
        else
        {
            rhs_.append(rhstemp_[i]);
        }
    }

    setThermo(thermoDatabase);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class GasThermoType, class SolidThermoType>
Foam::autoPtr<Foam::SurfaceReaction<GasThermoType, SolidThermoType>>
Foam::SurfaceReaction<GasThermoType, SolidThermoType>::New
(
    const speciesTable& gasSpecies,
    const speciesTable& species,    
    const HashPtrTable<GasThermoType>& gasThermoDatabase,
    const HashPtrTable<SolidThermoType>& thermoDatabase,    
    const dictionary& dict
)
{
    const word& reactionTypeName = dict.lookup("type");

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(reactionTypeName);

    // Backwards compatibility check. SurfaceReaction names used to 
    // have "SurfaceReaction" 
    // (SurfaceReaction<GasThermoType, SolidThermoType>::typeName_()) appended. 
    // This was removed as it is unnecessary given the context in which 
    // the reaction is specified. If this reaction name was not found, 
    // search also for the old name.
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        cstrIter = dictionaryConstructorTablePtr_->find
        (
            reactionTypeName.removeTrailing(typeName_())
        );
    }

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown reaction type "
            << reactionTypeName << nl << nl
            << "Valid reaction types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<SurfaceReaction<GasThermoType, SolidThermoType>>
    (
        cstrIter()(gasSpecies, species, gasThermoDatabase, thermoDatabase, dict)
    );
}


template<class GasThermoType, class SolidThermoType>
Foam::autoPtr<Foam::SurfaceReaction<GasThermoType, SolidThermoType>>
Foam::SurfaceReaction<GasThermoType, SolidThermoType>::New
(
    const speciesTable& gasSpecies,
    const speciesTable& species,    
    const HashPtrTable<GasThermoType>& gasThermoDatabase,
    const HashPtrTable<SolidThermoType>& thermoDatabase,    
    const objectRegistry& ob,
    const dictionary& dict
)
{
    // If the objectRegistry constructor table is empty
    // use the dictionary constructor table only
    if (!objectRegistryConstructorTablePtr_)
    {
        return New
        (
            gasSpecies, 
            species, 
            gasThermoDatabase, 
            thermoDatabase, 
            dict
        );
    }

    const word& reactionTypeName = dict.lookup("type");

    typename objectRegistryConstructorTable::iterator cstrIter =
        objectRegistryConstructorTablePtr_->find(reactionTypeName);

    // Backwards compatibility check. See above.
    if (cstrIter == objectRegistryConstructorTablePtr_->end())
    {
        cstrIter = objectRegistryConstructorTablePtr_->find
        (
            reactionTypeName.removeTrailing(typeName_())
        );
    }

    if (cstrIter == objectRegistryConstructorTablePtr_->end())
    {
        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(reactionTypeName);

        // Backwards compatibility check. See above.
        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            cstrIter = dictionaryConstructorTablePtr_->find
            (
                reactionTypeName.removeTrailing(typeName_())
            );
        }

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown reaction type "
                << reactionTypeName << nl << nl
                << "Valid reaction types are :" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << objectRegistryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<SurfaceReaction<GasThermoType, SolidThermoType>>
        (
            cstrIter()(gasSpecies, species, gasThermoDatabase, thermoDatabase, dict)
        );
    }

    return autoPtr<SurfaceReaction<GasThermoType, SolidThermoType>>
    (
        cstrIter()(gasSpecies, species, gasThermoDatabase, thermoDatabase, ob, dict)
    );
}


template<class GasThermoType, class SolidThermoType>
Foam::autoPtr<Foam::SurfaceReaction<GasThermoType, SolidThermoType>>
Foam::SurfaceReaction<GasThermoType, SolidThermoType>::New
(
    const speciesTable& gasSpecies,
    const speciesTable& species,    
    const PtrList<GasThermoType>& gasSpeciesThermo,
    const PtrList<SolidThermoType>& speciesThermo,    
    const dictionary& dict
)
{
    HashPtrTable<SolidThermoType> thermoDatabase;
    forAll(speciesThermo, i)
    {
        thermoDatabase.insert
        (
            speciesThermo[i].name(),
            speciesThermo[i].clone().ptr()
        );
    }

    HashPtrTable<GasThermoType> gasThermoDatabase;
    forAll(gasSpeciesThermo, i)
    {
        gasThermoDatabase.insert
        (
            gasSpeciesThermo[i].name(),
            gasSpeciesThermo[i].clone().ptr()
        );
    }

    return New(gasSpecies, species, gasThermoDatabase, thermoDatabase, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GasThermoType, class SolidThermoType>
void Foam::SurfaceReaction<GasThermoType, SolidThermoType>::write(Ostream& os) const
{
    OStringStream reaction;
    writeEntry
    (
        os,
        "reaction",
        catSpecieCoeffs::reactionStr
        (
            reaction, gasSpecies_, species_, lhsG_, rhsG_, lhs_, rhs_
        )        
    );
}


template<class GasThermoType, class SolidThermoType>
void Foam::SurfaceReaction<GasThermoType, SolidThermoType>::ddot
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& d
) const
{
}


template<class GasThermoType, class SolidThermoType>
void Foam::SurfaceReaction<GasThermoType, SolidThermoType>::fdot
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& f
) const
{
}


template<class GasThermoType, class SolidThermoType>
void Foam::SurfaceReaction<GasThermoType, SolidThermoType>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcdt,
    List<scalar> siteDensityList,
    List<label> siteNumberList    
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    scalar omegaI = omega
    (
        p, T, c, li, pf, cf, lRef, pr, cr, rRef, siteDensityList, siteNumberList
    );

    forAll(lhsG_, i)
    {
        const label gi  = lhsG_[i].index;
        const scalar gl = lhsG_[i].stoichCoeff;
        dcdt[gi] -= gl*omegaI; 
    }
    forAll(lhs_, i)
    {
        const label si  = lhs_[i].index + nG_;
        const scalar sl = lhs_[i].stoichCoeff;
        dcdt[si] -= sl*omegaI;
    }

    forAll(rhsG_, i)
    {
        const label gi  = rhsG_[i].index;
        const scalar gr = rhsG_[i].stoichCoeff;
        dcdt[gi] += gr*omegaI;
    }

    forAll(rhs_, i)
    {
        const label si = rhs_[i].index + nG_;
        const scalar sr = rhs_[i].stoichCoeff;
        dcdt[si] += sr*omegaI;
    }
}


template<class GasThermoType, class SolidThermoType>
Foam::scalar Foam::SurfaceReaction<GasThermoType, SolidThermoType>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef,
    List<scalar> siteDensityList,
    List<label> siteNumberList    
) const
{
    scalar clippedT = min(max(T, this->Tlow()), this->Thigh());

    //const scalar kf = this->kf(p, clippedT, c, li);
    const scalar kf = this->kf(p, clippedT, c, li, siteDensityList, siteNumberList);
    const scalar kr = this->kr(kf, p, clippedT, c, li);

    pf = 1;
    pr = 1;

    const label Nl = lhs_.size();
    const label Nr = rhs_.size();

    const label gNl = lhsG_.size(); //
    const label gNr = rhsG_.size(); //

    label slRef = 0;
    //lRef = lhs_[slRef].index;
    lRef = lhs_[slRef].index + nG_; //

    pf = kf;
    for (label s = 1; s < Nl; s++)
    {
        //const label si = lhs_[s].index;
        const label si = lhs_[s].index + nG_; //

        if (c[si] < c[lRef])
        {
            const scalar exp = lhs_[slRef].exponent;
            pf *= pow(max(c[lRef], 0), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = lhs_[s].exponent;
            pf *= pow(max(c[si], 0), exp);
        }
    }
    
    //
    for (label g = 0; g < gNl; g++)
    {
        const label gi = lhsG_[g].index;
        const scalar exp = lhsG_[g].exponent;
        pf *= pow(max(c[gi], 0), exp);
    }
    //

    cf = max(c[lRef], 0);

    {
        const scalar exp = lhs_[slRef].exponent;
        if (exp < 1)
        {
            if (cf > small)
            {
                pf *= pow(cf, exp - 1);
            }
            else
            {
                pf = 0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1);
        }
    }

    label srRef = 0;
    //rRef = rhs_[srRef].index;
    rRef = rhs_[srRef].index + nG_; //

    // Find the matrix element and element position for the rhs
    pr = kr;
    for (label s = 1; s < Nr; s++)
    {
        //const label si = rhs_[s].index;
        const label si = rhs_[s].index + nG_; //
        if (c[si] < c[rRef])
        {
            const scalar exp = rhs_[srRef].exponent;
            pr *= pow(max(c[rRef], 0), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = rhs_[s].exponent;
            pr *= pow(max(c[si], 0), exp);
        }
    }
 
   //
    for (label g = 0; g < gNr; g++)
    {
        const label gi = rhsG_[g].index;
        const scalar exp = rhsG_[g].exponent;
        pr *= pow(max(c[gi], 0), exp);
    }
   //
    cr = max(c[rRef], 0);

    {
        const scalar exp = rhs_[srRef].exponent;
        if (exp < 1)
        {
            if (cr > small)
            {
                pr *= pow(cr, exp - 1);
            }
            else
            {
                pr = 0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1);
        }
    }

    return pf*cf - pr*cr;
}


template<class GasThermoType, class SolidThermoType>
void Foam::SurfaceReaction<GasThermoType, SolidThermoType>::dwdc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarSquareMatrix& J,
    scalarField& dcdt,
    scalar& omegaI,
    scalar& kfwd,
    scalar& kbwd,
    const bool reduced,
    const List<label>& c2s,
    List<scalar> siteDensityList, 
    List<label> siteNumberList    
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    //omegaI = omega(p, T, c, li, pf, cf, lRef, pr, cr, rRef);
    omegaI = omega(p, T, c, li, pf, cf, lRef, pr, cr, rRef, siteDensityList, siteNumberList);

    forAll(lhs_, s)
    {
        label si = lhs_[s].index + nG_;
        const scalar sl = lhs_[s].stoichCoeff;
        dcdt[si] -= sl*omegaI;
    }
    forAll(lhsG_, g)
    {
        label gi = lhsG_[g].index;
        const scalar gl = lhsG_[g].stoichCoeff;
        dcdt[gi] -= gl*omegaI;
    }
    
    forAll(rhs_, s)
    {
        label si = rhs_[s].index + nG_;
        const scalar sr = rhs_[s].stoichCoeff;
        dcdt[si] += sr*omegaI;
    }
    forAll(rhsG_, g)
    {
        label gi = rhsG_[g].index;
        const scalar gr = rhsG_[g].stoichCoeff;
        dcdt[gi] += gr*omegaI;
    }

    //kfwd = this->kf(p, T, c, li);
    kfwd = this->kf(p, T, c, li, siteDensityList, siteNumberList);
    kbwd = this->kr(kfwd, p, T, c, li);

    forAll(lhs_, j)
    {
        const label sj = lhs_[j].index + nG_;
        scalar kf = kfwd;
        forAll(lhs_, i)
        {
            const label si = lhs_[i].index + nG_;
            const scalar el = lhs_[i].exponent;
            if (i == j)
            {
                if (el < 1)
                {
                    if (c[si] > SMALL)
                    {
                        kf *= el*pow(c[si] + VSMALL, el - 1);
                    }
                    else
                    {
                        kf = 0;
                    }
                }
                else
                {
                    kf *= el*pow(c[si], el - 1);
                }
            }
            else
            {
                kf *= pow(c[si], el);
            }
        }
        
        forAll(lhsG_, i)
        {
            const label si = lhsG_[i].index;
            const scalar el = lhsG_[i].exponent;
            kf *= pow(c[si], el);
        }

        forAll(lhs_, i)
        {
            const label si = lhs_[i].index + nG_;
            const scalar sl = lhs_[i].stoichCoeff;
            J(si, sj) -= sl*kf;
        }
        
        forAll(lhsG_, i)
        {
            const label si = lhsG_[i].index;
            const scalar sl = lhsG_[i].stoichCoeff;
            J(si, sj) -= sl*kf;
        }
                
        forAll(rhs_, i)
        {
            const label si = rhs_[i].index + nG_;
            const scalar sr = rhs_[i].stoichCoeff;
            J(si, sj) += sr*kf;
        }
        forAll(rhsG_, i)
        {
            const label si = rhsG_[i].index;
            const scalar sr = rhsG_[i].stoichCoeff;
            J(si, sj) += sr*kf;
        }
    }

    forAll(lhsG_, j)
    {
        const label sj = lhsG_[j].index;
        scalar kf = kfwd;
        forAll(lhsG_, i)
        {
            const label si = lhsG_[i].index;
            const scalar el = lhsG_[i].exponent;
            if (i == j)
            {
                if (el < 1)
                {
                    if (c[si] > SMALL)
                    {
                        kf *= el*pow(c[si] + VSMALL, el - 1);
                    }
                    else
                    {
                        kf = 0;
                    }
                }
                else
                {
                    kf *= el*pow(c[si], el - 1);
                }
            }
            else
            {
                kf *= pow(c[si], el);
            }
        }
        
        forAll(lhs_, i)
        {
            const label si = lhs_[i].index + nG_;
            const scalar el = lhs_[i].exponent;
            kf *= pow(c[si], el);
        }

        forAll(lhs_, i)
        {
            const label si = lhs_[i].index + nG_;
            const scalar sl = lhs_[i].stoichCoeff;
            J(si, sj) -= sl*kf;
        }
        
        forAll(lhsG_, i)
        {
            const label si = lhsG_[i].index;
            const scalar sl = lhsG_[i].stoichCoeff;
            J(si, sj) -= sl*kf;
        }
        
        forAll(rhs_, i)
        {
            const label si = rhs_[i].index + nG_;
            const scalar sr = rhs_[i].stoichCoeff;
            J(si, sj) += sr*kf;
        }
        forAll(rhsG_, i)
        {
            const label si = rhsG_[i].index;
            const scalar sr = rhsG_[i].stoichCoeff;
            J(si, sj) += sr*kf;
        }
    }


    forAll(rhs_, j)
    {
        const label sj = rhs_[j].index + nG_;
        scalar kr = kbwd;
        forAll(rhs_, i)
        {
            const label si = rhs_[i].index + nG_;
            const scalar er = rhs_[i].exponent;
            if (i == j)
            {
                if (er < 1)
                {
                    if (c[si] > SMALL)
                    {
                        kr *= er*pow(c[si] + VSMALL, er - 1);
                    }
                    else
                    {
                        kr = 0;
                    }
                }
                else
                {
                    kr *= er*pow(c[si], er - 1);
                }
            }
            else
            {
                kr *= pow(c[si], er);
            }
        }
        
        forAll(rhsG_, i)
        {
            const label si = rhsG_[i].index;
            const scalar er = rhsG_[i].exponent;
            kr *= pow(c[si], er);
        }

        forAll(lhs_, i)
        {
            const label si = lhs_[i].index + nG_;
            const scalar sl = lhs_[i].stoichCoeff;
            J(si, sj) += sl*kr;
        }
        
        forAll(lhsG_, i)
        {
            const label si = lhsG_[i].index;
            const scalar sl = lhsG_[i].stoichCoeff;
            J(si, sj) += sl*kr;
        }
        
        forAll(rhs_, i)
        {
            const label si = rhs_[i].index + nG_;
            const scalar sr = rhs_[i].stoichCoeff;
            J(si, sj) -= sr*kr;
        }
        
        forAll(rhsG_, i)
        {
            const label si = rhsG_[i].index;
            const scalar sr = rhsG_[i].stoichCoeff;
            J(si, sj) -= sr*kr;
        }
    }

    forAll(rhsG_, j)
    {
        const label sj = rhsG_[j].index;
        scalar kr = kbwd;
        forAll(rhsG_, i)
        {
            const label si = rhsG_[i].index;
            const scalar er = rhsG_[i].exponent;
            if (i == j)
            {
                if (er < 1)
                {
                    if (c[si] > SMALL)
                    {
                        kr *= er*pow(c[si] + VSMALL, er - 1);
                    }
                    else
                    {
                        kr = 0;
                    }
                }
                else
                {
                    kr *= er*pow(c[si], er - 1);
                }
            }
            else
            {
                kr *= pow(c[si], er);
            }
        }
        
        forAll(rhs_, i)
        {
            const label si = rhs_[i].index + nG_;
            const scalar er = rhs_[i].exponent;
            kr *= pow(c[si], er);
        }

        forAll(lhs_, i)
        {
            const label si = lhs_[i].index + nG_;
            const scalar sl = lhs_[i].stoichCoeff;
            J(si, sj) += sl*kr;
        }
        
        forAll(lhsG_, i)
        {
            const label si = lhsG_[i].index;
            const scalar sl = lhsG_[i].stoichCoeff;
            J(si, sj) += sl*kr;
        }
        
        forAll(rhs_, i)
        {
            const label si = rhs_[i].index + nG_;
            const scalar sr = rhs_[i].stoichCoeff;
            J(si, sj) -= sr*kr;
        }
        forAll(rhsG_, i)
        {
            const label si = rhsG_[i].index;
            const scalar sr = rhsG_[i].stoichCoeff;
            J(si, sj) -= sr*kr;
        }
    }

}


template<class GasThermoType, class SolidThermoType>
void Foam::SurfaceReaction<GasThermoType, SolidThermoType>::dwdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const scalar omegaI,
    const scalar kfwd,
    const scalar kbwd,
    scalarSquareMatrix& J,
    const bool reduced,
    const List<label>& c2s,
    const label indexT
) const
{
    scalar kf = kfwd;
    scalar kr = kbwd;

    scalar dkfdT = this->dkfdT(p, T, c, li);
    scalar dkrdT = this->dkrdT(p, T, c, li, dkfdT, kr);

    scalar sumExp = 0.0;
    forAll(lhs_, i)
    {
        const label si = lhs_[i].index;
        const scalar el = lhs_[i].exponent;
        const scalar cExp = pow(c[si], el);
        dkfdT *= cExp;
        kf *= cExp;
        sumExp += el;
    }
    kf *= -sumExp/T;

    sumExp = 0.0;
    forAll(rhs_, i)
    {
        const label si = rhs_[i].index;
        const scalar er = rhs_[i].exponent;
        const scalar cExp = pow(c[si], er);
        dkrdT *= cExp;
        kr *= cExp;
        sumExp += er;
    }
    kr *= -sumExp/T;

    // dqidT includes the third-body (or pressure dependent) effect
    scalar dqidT = dkfdT - dkrdT + kf - kr;

    // For reactions including third-body efficiencies or pressure dependent
    // reaction, an additional term is needed
    scalar dcidT = this->dcidT(p, T, c, li);
    dcidT *= omegaI;

    // J(i, indexT) = sum_reactions nu_i dqdT
    forAll(lhs_, i)
    {
        const label si = reduced ? c2s[lhs_[i].index] : lhs_[i].index;
        const scalar sl = lhs_[i].stoichCoeff;
        J(si, indexT) -= sl*(dqidT + dcidT);
    }
    forAll(rhs_, i)
    {
        const label si = reduced ? c2s[rhs_[i].index] : rhs_[i].index;
        const scalar sr = rhs_[i].stoichCoeff;
        J(si, indexT) += sr*(dqidT + dcidT);
    }
}


template<class GasThermoType, class SolidThermoType>
const Foam::speciesTable& Foam::SurfaceReaction<GasThermoType, SolidThermoType>::
species() const
{
    return species_;
}


template<class GasThermoType, class SolidThermoType>
const Foam::speciesTable& Foam::SurfaceReaction<GasThermoType, SolidThermoType>::
gasSpecies() const
{
    return gasSpecies_;
}


// ************************************************************************* //
