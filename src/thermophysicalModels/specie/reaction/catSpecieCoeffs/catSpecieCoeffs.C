/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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

#include "catSpecieCoeffs.H"
#include "DynamicList.H"
#include "OStringStream.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::catSpecieCoeffs::catSpecieCoeffs
(
    const speciesTable& species,
    const speciesTable& solidSpecies, //for cat    
    Istream& is
)
{
    token t(is);
    if (t.isNumber())
    {
        stoichCoeff = t.number();
        is >> t;
    }
    else
    {
        stoichCoeff = 1;
    }

    exponent = stoichCoeff;

    if (t.isWord())
    {
        word specieName = t.wordToken();

        size_t i = specieName.find('^');

        if (i != word::npos)
        {
            string exponentStr = specieName
            (
                i + 1,
                specieName.size() - i - 1
            );
            exponent = atof(exponentStr.c_str());
            specieName = specieName(0, i);
        }

        if (species.found(specieName))
        {
            index = species[specieName];
            isgas = true; //for cat            
        }
        else if (solidSpecies.found(specieName)) //for cat
        {
            index = solidSpecies[specieName];
            isgas = false; 
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "Specie " << specieName
                << " not found in table " << species
                << exit(FatalIOError);

            index = -1;
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Expected a word but found " << t.info()
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::catSpecieCoeffs::setLRhs
(
    Istream& is,
    const speciesTable& species,
    List<catSpecieCoeffs>& lhs,
    List<catSpecieCoeffs>& rhs
)
{
    DynamicList<catSpecieCoeffs> dlrhs;

    while (is.good())
    {
        //dlrhs.append(catSpecieCoeffs(species, is));
        dlrhs.append(catSpecieCoeffs(species, species, is)); //for cat
        if (dlrhs.last().index != -1)
        {
            token t(is);
            if (t.isPunctuation())
            {
                if (t == token::ADD)
                {
                }
                else if (t == token::ASSIGN)
                {
                    lhs = dlrhs.shrink();
                    dlrhs.clear();
                }
                else
                {
                    rhs = dlrhs.shrink();
                    is.putBack(t);
                    return;
                }
            }
            else
            {
                rhs = dlrhs.shrink();
                is.putBack(t);
                return;
            }
        }
        else
        {
            dlrhs.remove();
            if (is.good())
            {
                token t(is);
                if (t.isPunctuation())
                {
                    if (t == token::ADD)
                    {
                    }
                    else if (t == token::ASSIGN)
                    {
                        lhs = dlrhs.shrink();
                        dlrhs.clear();
                    }
                    else
                    {
                        rhs = dlrhs.shrink();
                        is.putBack(t);
                        return;
                    }
                }
            }
            else
            {
                if (!dlrhs.empty())
                {
                    rhs = dlrhs.shrink();
                }
                return;
            }
        }
    }

    FatalIOErrorInFunction(is)
        << "Cannot continue reading reaction data from stream"
        << exit(FatalIOError);
}


void Foam::catSpecieCoeffs::reactionStr
(
    OStringStream& reaction,
    const speciesTable& species,
    const List<catSpecieCoeffs>& scs
)
{
    for (label i = 0; i < scs.size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(scs[i].stoichCoeff - 1) > small)
        {
            reaction << scs[i].stoichCoeff;
        }
        reaction << species[scs[i].index];
        if (mag(scs[i].exponent - scs[i].stoichCoeff) > small)
        {
            reaction << "^" << scs[i].exponent;
        }
    }
}


Foam::string Foam::catSpecieCoeffs::reactionStr
(
    OStringStream& reaction,
    const speciesTable& species,
    const List<catSpecieCoeffs>& lhs,
    const List<catSpecieCoeffs>& rhs
)
{
    reactionStr(reaction, species, lhs);
    reaction << " = ";
    reactionStr(reaction, species, rhs);
    return reaction.str();
}


// for cat
void Foam::catSpecieCoeffs::setLRhsCat
(
    Istream& is,
    const speciesTable& species,
    const speciesTable& solidSpecies,
    List<catSpecieCoeffs>& lhs,
    List<catSpecieCoeffs>& rhs
)
{
    DynamicList<catSpecieCoeffs> dlrhs;

    while (is.good())
    {
        //dlrhs.append(catSpecieCoeffs(species, is));
        dlrhs.append(catSpecieCoeffs(species, solidSpecies, is)); //for cat

        if (dlrhs.last().index != -1)
        {
            token t(is);
            if (t.isPunctuation())
            {
                if (t == token::ADD)
                {
                }
                else if (t == token::ASSIGN)
                {
                    lhs = dlrhs.shrink();
                    dlrhs.clear();
                }
                else
                {
                    rhs = dlrhs.shrink();
                    is.putBack(t);
                    return;
                }
            }
            else
            {
                rhs = dlrhs.shrink();
                is.putBack(t);
                return;
            }
        }
        else
        {
            dlrhs.remove();
            if (is.good())
            {
                token t(is);
                if (t.isPunctuation())
                {
                    if (t == token::ADD)
                    {
                    }
                    else if (t == token::ASSIGN)
                    {
                        lhs = dlrhs.shrink();
                        dlrhs.clear();
                    }
                    else
                    {
                        rhs = dlrhs.shrink();
                        is.putBack(t);
                        return;
                    }
                }
            }
            else
            {
                if (!dlrhs.empty())
                {
                    rhs = dlrhs.shrink();
                }
                return;
            }
        }
    }

    FatalIOErrorInFunction(is)
        << "Cannot continue reading reaction data from stream"
        << exit(FatalIOError);
}


void Foam::catSpecieCoeffs::reactionStr
(
    OStringStream& reaction,
    const speciesTable& gasSpecies,
    const speciesTable& species,    
    const List<catSpecieCoeffs>& scsG,
    const List<catSpecieCoeffs>& scs    
)
{
    for (label i = 0; i < scsG.size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(scsG[i].stoichCoeff - 1) > small)
        {
            reaction << scsG[i].stoichCoeff;
        }
        reaction << gasSpecies[scsG[i].index];
        if (mag(scsG[i].exponent - scsG[i].stoichCoeff) > small)
        {
            reaction << "^" << scsG[i].exponent;
        }
    }

    for (label i = 0; i < scs.size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(scs[i].stoichCoeff - 1) > small)
        {
            reaction << scs[i].stoichCoeff;
        }
        reaction << species[scs[i].index];
        if (mag(scs[i].exponent - scs[i].stoichCoeff) > small)
        {
            reaction << "^" << scs[i].exponent;
        }
    }
}


Foam::string Foam::catSpecieCoeffs::reactionStr
(
    OStringStream& reaction,
    const speciesTable& gasSpecies,
    const speciesTable& species,    
    const List<catSpecieCoeffs>& lhsG,
    const List<catSpecieCoeffs>& rhsG,
    const List<catSpecieCoeffs>& lhs,
    const List<catSpecieCoeffs>& rhs    
)
{
    reactionStr(reaction, gasSpecies, species, lhsG, lhs);
    reaction << " = ";
    reactionStr(reaction, gasSpecies, species, rhsG, rhs);
    return reaction.str();
}

// 

// ************************************************************************* //
