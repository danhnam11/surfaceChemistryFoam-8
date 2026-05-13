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

Class
    Foam::CoefficientsManager

Description
    Transport package using Standard Kinetic Theory model. The mixture is 
    calculated based on the modification of the Wilke semi-empirical formulas. 

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Mass Diffusivity is based on mixture averaged model in which binary diffusion 
    coefficients are obtained from Standard Kinetic Theory model.

    It has been validated against NIST and real-fluid based OPPDIF. 
 
    For details of equations in the model, refer to TRANSPORT program 
    in CHEMKIN II @ R. J. Kee et al. CHEMKIN collection, Release 3.6, 
    Reaction Design, Inc., San Diego, CA (2001). 

SourceFiles
    CoefficientsManager.H
    CoefficientsManager.C

\*---------------------------------------------------------------------------*/

#include "CoefficientsManager.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::CoefficientsManager::CoefficientsManager()
:
    muCoeffs_(),
    kappaCoeffs_(),
    DijCoeffs_()
{
}

Foam::CoefficientsManager::CoefficientsManager(const Foam::dictionary& dict)
:
    muCoeffs_(),
    kappaCoeffs_(),
    DijCoeffs_()
{
    dict.lookup("muCoeffs") >> muCoeffs_;
    dict.lookup("kappaCoeffs") >> kappaCoeffs_;
    dict.lookup("DijCoeffs") >> DijCoeffs_;

    label nSpecies = DijCoeffs_.size();
    printCoeffs();

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CoefficientsManager::setMkCoeffs
(
    const List<List<scalar>>& muCoeffsMk,
    const List<List<scalar>>& kappaCoeffsMk,
    const List<List<List<scalar>>>& DijCoeffsMk
)
{
    muCoeffsMk_ = muCoeffsMk;
    kappaCoeffsMk_ = kappaCoeffsMk;
    DijCoeffsMk_ = DijCoeffsMk;
}

void Foam::CoefficientsManager::setMuKappaMkCoeffs
(
    const List<List<scalar>>& muCoeffsMk,
    const List<List<scalar>>& kappaCoeffsMk
)
{
    muCoeffsMk_ = muCoeffsMk;
    kappaCoeffsMk_ = kappaCoeffsMk;
}


void Foam::CoefficientsManager::printCoeffs() const
{
    Info << "muCoeffs: " << muCoeffs_ << nl;
    Info << "kappaCoeffs: " << kappaCoeffs_ << nl;
    Info << "DijCoeffs: " << DijCoeffs_ << nl;
}

const Foam::List<Foam::scalar>& Foam::CoefficientsManager::muCoeffs() const
{
    return muCoeffs_; 
}

const Foam::List<Foam::scalar>& Foam::CoefficientsManager::kappaCoeffs() const
{ 
    return kappaCoeffs_; 
}

const Foam::List<Foam::List<Foam::scalar>>& Foam::CoefficientsManager::DijCoeffs() const
{ 
    return DijCoeffs_; 
}

const Foam::List<Foam::List<Foam::scalar>>& Foam::CoefficientsManager::muCoeffsMk() const
{
    return muCoeffsMk_;
}

const Foam::List<Foam::List<Foam::scalar>>& Foam::CoefficientsManager::kappaCoeffsMk() const
{
    return kappaCoeffsMk_;
}

const Foam::List<Foam::List<Foam::List<Foam::scalar>>>& Foam::CoefficientsManager::DijCoeffsMk() const
{
    return DijCoeffsMk_;
}


// ************************************************************************* //
