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

inline Foam::StickArrheniusSurfaceReactionRate::StickArrheniusSurfaceReactionRate
(
    const scalar a,
    const scalar b,
    const scalar c
)
:
    //A_(A),
    //beta_(beta),
    //Ta_(Ta)
    a_(a),
    b_(b),
    c_(c)    
{}


inline Foam::StickArrheniusSurfaceReactionRate::StickArrheniusSurfaceReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    //A_(dict.lookup<scalar>("A")),
    //beta_(dict.lookup<scalar>("beta")),
    //Ta_(dict.lookup<scalar>("Ta"))
    a_(dict.lookup<scalar>("A")),
    b_(dict.lookup<scalar>("beta")),
    c_(dict.lookup<scalar>("Ta"))    
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::StickArrheniusSurfaceReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    scalar ak = a_;

    if (mag(b_) > vSmall)
    {
        ak *= pow(T, b_);
    }

    if (mag(c_) > vSmall)
    {
        ak *= exp(-c_/T);
    }

    return min(1, ak);
}


inline Foam::scalar Foam::StickArrheniusSurfaceReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField&,
    const label
) const
{
    scalar ak = a_;

    if (mag(b_) > vSmall)
    {
        ak *= pow(T, b_);
    }

    if (mag(c_) > vSmall)
    {
        ak *= exp(-c_/T);
    }

    return ak*(b_+c_/T)/T;
}

/*
inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalar>>&
Foam::StickArrheniusSurfaceReactionRate::beta() const
{
    return NullObjectRef<List<Tuple2<label, scalar>>>();
}
*/

inline const Foam::List<Foam::Tuple2<Foam::label, Foam::scalarField>>&
Foam::StickArrheniusSurfaceReactionRate::coverageModificationCoeffs() const
{
    return NullObjectRef<List<Tuple2<label, scalarField>>>();
}



inline void Foam::StickArrheniusSurfaceReactionRate::dcidc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcidc
) const
{}


inline Foam::scalar Foam::StickArrheniusSurfaceReactionRate::dcidT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return 0;
}


inline void Foam::StickArrheniusSurfaceReactionRate::write(Ostream& os) const
{
    writeEntry(os, "A", a_);
    writeEntry(os, "beta", b_);
    writeEntry(os, "Ta", c_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const StickArrheniusSurfaceReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
