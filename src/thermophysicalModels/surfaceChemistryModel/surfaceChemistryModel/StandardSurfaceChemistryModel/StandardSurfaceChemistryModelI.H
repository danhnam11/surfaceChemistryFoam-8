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

#include "volFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
inline Foam::label
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
nEqns() const
{
    // nEqns = number of species + temperature + pressure
    return nSpecie_ + 2;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
RR()
{
    return RRs_;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
inline const Foam::PtrList<Foam::SurfaceReaction<GasThermoType, SolidThermoType>>&
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
reactions() const
{
    return reactions_;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
inline const Foam::PtrList<GasThermoType>&
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
gasSpecieThermos() const
{
    return gasSpecieThermos_;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
inline const Foam::PtrList<SolidThermoType>&
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
specieThermos() const
{
    return specieThermos_;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
inline Foam::label
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
nSpecie() const
{
    return nSpecie_;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
inline Foam::label
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
nGasSpecie() const
{
    return nGasSpecie_;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
inline Foam::label
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
nSolidSpecie() const
{
    return nSolidSpecie_;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
inline Foam::label
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
nReaction() const
{
    return nReaction_;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
inline Foam::scalar
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
Treact() const
{
    return Treact_;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
inline Foam::scalar&
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
Treact()
{
    return Treact_;
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
RR
(
    const label i
) const
{
    return RRs_[i];
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
RR
(
    const label i
)
{
    return RRs_[i];
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
RRg
(
    const label i
) const
{
    return RRg_[i];
}


template
<
    class GasReactionThermo, 
    class SolidReactionThermo,     
    class GasThermoType,
    class SolidThermoType    
>
Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::StandardSurfaceChemistryModel
<GasReactionThermo, SolidReactionThermo, GasThermoType, SolidThermoType>::
RRg
(
    const label i
)
{
    return RRg_[i];
}



// ************************************************************************* //
