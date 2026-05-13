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

inline Foam::ArrheniusSurfaceReactionRate::ArrheniusSurfaceReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta
)
:
    A_(A),
    beta_(beta),
    Ta_(Ta)
{}


inline Foam::ArrheniusSurfaceReactionRate::ArrheniusSurfaceReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    A_(dict.lookup<scalar>("A")),
    beta_(dict.lookup<scalar>("beta")),
    Ta_(dict.lookup<scalar>("Ta"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::ArrheniusSurfaceReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    scalar ak = A_;

    if (mag(beta_) > vSmall)
    {
        ak *= pow(T, beta_);
    }

    if (mag(Ta_) > vSmall)
    {
        ak *= exp(-Ta_/T);
    }

    return ak;
}


inline Foam::scalar Foam::ArrheniusSurfaceReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    scalar ak = A_;

    if (mag(beta_) > vSmall)
    {
        ak *= pow(T, beta_);
    }

    if (mag(Ta_) > vSmall)
    {
        ak *= exp(-Ta_/T);
    }

    return ak*(beta_+Ta_/T)/T;
}

/*
inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::ArrheniusSurfaceReactionRate::beta() const
{
    return NullObjectRef<List<Tuple2<label, scalar>>>();
}
*/

inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalarField>>&
Foam::ArrheniusSurfaceReactionRate::coverageModificationCoeffs() const
{
    return NullObjectRef<List<Tuple2<label, scalarField>>>();
}



inline void Foam::ArrheniusSurfaceReactionRate::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcidc
) const
{}


inline Foam::scalar Foam::ArrheniusSurfaceReactionRate::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return 0;
}


inline void Foam::ArrheniusSurfaceReactionRate::write(Ostream& os) const
{
    writeEntry(os, "A", A_);
    writeEntry(os, "beta", beta_);
    writeEntry(os, "Ta", Ta_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ArrheniusSurfaceReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
