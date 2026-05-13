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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::coverageModificationArrheniusSurfaceReactionRate::
coverageModificationArrheniusSurfaceReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const coverageModifications& covM
)
:
    ArrheniusSurfaceReactionRate(A, beta, Ta),
    coverageModifications_(covM),
    coverageModificationCoeffs_(covM.covCoeffs())    
{}


inline Foam::coverageModificationArrheniusSurfaceReactionRate::
coverageModificationArrheniusSurfaceReactionRate
(
    const speciesTable& species,
    const dictionary& dict
)
:
    ArrheniusSurfaceReactionRate
    (
        species,
        dict
    ),
    coverageModifications_(species, dict),
    coverageModificationCoeffs_(coverageModifications_.covCoeffs())    
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::coverageModificationArrheniusSurfaceReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return ArrheniusSurfaceReactionRate::operator()(p, T, c, li);
}


inline Foam::scalar Foam::coverageModificationArrheniusSurfaceReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return ArrheniusSurfaceReactionRate::ddT(p, T, c, li);
}


inline void Foam::coverageModificationArrheniusSurfaceReactionRate::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcidc
) const
{}


inline Foam::scalar Foam::coverageModificationArrheniusSurfaceReactionRate::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return -1.0/T;
}


inline void Foam::coverageModificationArrheniusSurfaceReactionRate::write(Ostream& os) const
{
    ArrheniusSurfaceReactionRate::write(os);
    coverageModifications_.write(os);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const coverageModificationArrheniusSurfaceReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
