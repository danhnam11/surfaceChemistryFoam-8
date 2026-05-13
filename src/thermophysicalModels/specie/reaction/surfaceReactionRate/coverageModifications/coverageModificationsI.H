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

#include "Tuple2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::coverageModifications::coverageModifications
(
    const speciesTable& species,
    const dictionary& dict
)
:
    species_(species)
{
    if (dict.found("coeffs"))
    {
        List<Tuple2<word, scalarField>> coeffs(dict.lookup("coeffs"));
        if (coeffs.size() != species_.size())
        {
            FatalIOErrorInFunction(dict)
                << "number of coverage modifications = " << coeffs.size()
                << " is not equal to the number of species " << species_.size()
                << exit(FatalIOError);
        }

        forAll(coeffs, i)
        {   
            covCoeffs_.append
            (   
                Tuple2<label, scalarField>
                (   
                    species_[coeffs[i].first()], 
                    coeffs[i].second()
                )
            );

        }
    }
}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //

inline void Foam::coverageModifications::write(Ostream& os) const
{
    List<Tuple2<word, scalarField>> coeffs(species_.size());
    //List<Tuple2<word, scalar>> coeffs(species_.size());
    forAll(coeffs, i)
    {
        coeffs[i].first() = species_[i];
        coeffs[i].second() = covCoeffs_[i].second();
    }

    writeEntry(os, "coeffs", coeffs);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const coverageModifications& covM
)
{
    covM.write(os);
    return os;
}


// ************************************************************************* //
