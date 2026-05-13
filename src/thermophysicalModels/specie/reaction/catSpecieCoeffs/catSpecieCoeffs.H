/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::catSpecieCoeffs

    developed by Danh Nam Nguyen
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).

Description
    Class to hold the specie index and its coefficients in the
    reaction rate expression for surface reactions. 

SourceFiles
    catSpecieCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef catSpecieCoeffs_H
#define catSpecieCoeffs_H

#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class catSpecieCoeffs;
Ostream& operator<<(Ostream&, const catSpecieCoeffs&);


/*---------------------------------------------------------------------------*\
                         Class catSpecieCoeffs Declaration
\*---------------------------------------------------------------------------*/

class catSpecieCoeffs
{
public:

    label index;
    scalar stoichCoeff;
    scalar exponent;
    bool isgas; //for cat    

    catSpecieCoeffs()
    :
        index(-1),
        stoichCoeff(0),
        exponent(1),
        isgas(true) //for cat        
    {}

    // catSpecieCoeffs(const speciesTable& species, Istream& is);
    catSpecieCoeffs
    (
        const speciesTable& species, 
        const speciesTable& solidSpecies,
        Istream& is
    );


    bool operator==(const catSpecieCoeffs& sc) const
    {
        return index == sc.index;
    }

    bool operator!=(const catSpecieCoeffs& sc) const
    {
        return index != sc.index;
    }

    friend Ostream& operator<<(Ostream& os, const catSpecieCoeffs& sc)
    {
        os  << sc.index << token::SPACE
            << sc.stoichCoeff << token::SPACE
            //<< sc.exponent;
            << sc.exponent << token::SPACE << sc.isgas; // for cat          
        return os;
    }

    //- Construct the left- and right-hand-side reaction coefficients
    static void setLRhs
    (
        Istream&,
        const speciesTable&,
        List<catSpecieCoeffs>& lhs,
        List<catSpecieCoeffs>& rhs
    );

    //- Construct the left- and right-hand-side reaction coefficients for catalyst
    static void setLRhsCat //for cat
    (
        Istream&,
        const speciesTable&,
        const speciesTable&,
        List<catSpecieCoeffs>& lhs,
        List<catSpecieCoeffs>& rhs
    );


    //- Write the string representation of the catSpecieCoeffs list
    static void reactionStr
    (
        OStringStream& reaction,
        const speciesTable&,
        const List<catSpecieCoeffs>& scs
    );

    //- Return string representation of reaction
    static string reactionStr
    (
        OStringStream& reaction,
        const speciesTable&,
        const List<catSpecieCoeffs>& lhs,
        const List<catSpecieCoeffs>& rhs
    );

    //- Write the string representation of the catSpecieCoeffs list
    static void reactionStr // for cat
    (
        OStringStream& reaction,
        const speciesTable&,
        const speciesTable&,        
        const List<catSpecieCoeffs>& scsG,
        const List<catSpecieCoeffs>& scs        
    );

    //- Return string representation of reaction
    static string reactionStr //for cat
    (
        OStringStream& reaction,
        const speciesTable&,
        const speciesTable&,
        const List<catSpecieCoeffs>& lhsG,
        const List<catSpecieCoeffs>& rhsG,                
        const List<catSpecieCoeffs>& lhs,
        const List<catSpecieCoeffs>& rhs
    );


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
