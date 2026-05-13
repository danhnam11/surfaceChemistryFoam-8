/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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

#include "catThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(catThermo, 0);
    defineRunTimeSelectionTable(catThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::catThermo::catThermo(const fvMesh& mesh, const word& phaseName)
:
    basicCatThermo(mesh, phaseName),
    p_(lookupOrConstruct(mesh, "p.catalyst"))
{}



Foam::catThermo::catThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    basicCatThermo(mesh, dict, phaseName),
    p_(lookupOrConstruct(mesh, "p.catalyst"))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::catThermo> Foam::catThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicCatThermo::New<catThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::catThermo::~catThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField& Foam::catThermo::p()
{
    return p_;
}


const Foam::volScalarField& Foam::catThermo::p() const
{
    return p_;
}


Foam::tmp<Foam::volScalarField> Foam::catThermo::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::scalarField> Foam::catThermo::nu(const label patchi) const
{
    return mu(patchi)/rho(patchi);
}


// ************************************************************************* //
