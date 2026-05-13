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

Class
    Foam::multiComponentCatMixture

    developed by Danh Nam Nguyen
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).

Description

    Foam::multiComponentCatMixture
    This is original multiComponentCatMixture class but has been modified 
    to be compatible with interface basicMultiComponentCatMixture 
    to avoid using of hardcoding in StandardChemistryModel class 

SourceFiles
    multiComponentCatMixture.C

\*---------------------------------------------------------------------------*/

#ifndef multiComponentCatMixture_H
#define multiComponentCatMixture_H

#include "basicMultiComponentCatMixture.H"
#include "HashPtrTable.H"
#include "specieElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiComponentCatMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class multiComponentCatMixture
:
    public basicMultiComponentCatMixture<ThermoType>
{
    // Private Data

        //- Temporary storage for the cell/face mixture thermo data
        mutable ThermoType mixture_;

        //- Temporary storage for the volume weighted
        //  cell/face mixture thermo data
        mutable ThermoType mixtureVol_;

    // Private Member Functions

        //- Correct the mass fractions to sum to 1
        void correctMassFractions();

public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        multiComponentCatMixture(const dictionary&, const fvMesh&, const word&);

        //- Disallow default bitwise copy construction
        multiComponentCatMixture
        (
            const multiComponentCatMixture<ThermoType>&
        ) = delete;


    //- Destructor
    virtual ~multiComponentCatMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "multiComponentCatMixture<" + ThermoType::typeName() + '>';
        }

        const ThermoType& cellMixture(const label celli) const;

        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        const ThermoType& cellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const;

        const ThermoType& patchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const;

        //- Read dictionary
        void read(const dictionary& thermoDict)
        {
            basicMultiComponentCatMixture<ThermoType>::read(thermoDict);
        }

        // - Return the number of species in the mixture - Nam
        inline const label& numberOfSpecies() const
        {
            return basicMultiComponentCatMixture<ThermoType>::numberOfSpecies();
        }
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multiComponentCatMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
