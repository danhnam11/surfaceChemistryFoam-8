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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GasThermoType, class SolidThermoType>
inline const word& SurfaceReaction<GasThermoType, SolidThermoType>::name() const
{
    return name_;
}


template<class GasThermoType, class SolidThermoType>
inline scalar SurfaceReaction<GasThermoType, SolidThermoType>::Tlow() const
{
    return Tlow_;
}


template<class GasThermoType, class SolidThermoType>
inline scalar SurfaceReaction<GasThermoType, SolidThermoType>::Thigh() const
{
    return Thigh_;
}


template<class GasThermoType, class SolidThermoType>
inline const List<catSpecieCoeffs>& SurfaceReaction<GasThermoType, SolidThermoType>::
lhs() const
{
    return lhs_;
}


template<class GasThermoType, class SolidThermoType>
inline const List<catSpecieCoeffs>& SurfaceReaction<GasThermoType, SolidThermoType>::
rhs() const
{
    return rhs_;
}


template<class GasThermoType, class SolidThermoType>
inline const List<catSpecieCoeffs>& SurfaceReaction<GasThermoType, SolidThermoType>::
lhsG() const
{
    return lhsG_;
}


template<class GasThermoType, class SolidThermoType>
inline const List<catSpecieCoeffs>& SurfaceReaction<GasThermoType, SolidThermoType>::
rhsG() const
{
    return rhsG_;
}



// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class GasThermoType, class SolidThermoType>
inline Ostream& operator<<(Ostream& os, const SurfaceReaction<GasThermoType, SolidThermoType>& r)
{
    r.write(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
