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
    Foam::DTMMultiComponentMixture

Description

    developed by:
    Nam Danh Nguyen and Jae Hun Lee
    Advisor: Prof. Chun Sang Yoo 
    Clean combustion & energy research (CCER) Lab - Dept. of Mech. Engineering
    Ulsan Institute of Science and Technology (UNIST) - Korea

SourceFiles
    DTMMultiComponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef DTMMultiComponentMixture_H
#define DTMMultiComponentMixture_H

#include "basicMultiComponentMixture.H"
#include "HashPtrTable.H"
#include "specieElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class DTMMultiComponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class DTMMultiComponentMixture
:
    public basicMultiComponentMixture<ThermoType>
{
    // Private Data

        //- Temporary storage for the cell/face mixture thermo data
        mutable ThermoType mixture_;

        //- Temporary storage for the volume weighted
        //  cell/face mixture thermo data
        mutable ThermoType mixtureVol_;

      //- For mass diffusivity in Standard Kinetic Theory model
        //- List of species molecular weight
        List<scalar> ListW_;

        List<scalar> linearityMk_;
        List<scalar> epsilonOverKbMk_;
        List<scalar> sigmaMk_;
        List<scalar> miuiMk_;
        List<scalar> polarMk_;
        List<scalar> ZrotMk_;
        List<List<List<scalar>>> CpCoeffTableMk_;

        List<List<scalar>> EPSILONijOVERKB_;
        List<List<scalar>> DELTAij_;
        List<List<scalar>> Mij_;
        List<List<scalar>> SIGMAij_;

        const fvMesh& mesh_;


    // Private Member Functions

        //- Correct the mass fractions to sum to 1
        void correctMassFractions();


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        DTMMultiComponentMixture(const dictionary&, const fvMesh&, const word&);

        //- Disallow default bitwise copy construction
        DTMMultiComponentMixture
        (
            const DTMMultiComponentMixture<ThermoType>&
        ) = delete;


    //- Destructor
    virtual ~DTMMultiComponentMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "DTMMultiComponentMixture<" + ThermoType::typeName() + '>';
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
            basicMultiComponentMixture<ThermoType>::read(thermoDict);
        }

        // - Return the number of species in the mixture - Nam
        inline const label& numberOfSpecies() const
        {
            return basicMultiComponentMixture<ThermoType>::numberOfSpecies();
        }
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DTMMultiComponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
