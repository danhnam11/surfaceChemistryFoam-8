/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2020 OpenFOAM Foundation
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
    Foam::basicSpecieCatMixture

    developed by Danh Nam Nguyen
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).

Description
    Specialization of basicCatMixture for a mixture consisting of a number for
    molecular species.

SourceFiles
    basicSpecieCatMixture.C

\*---------------------------------------------------------------------------*/

#ifndef basicSpecieCatMixture_H
#define basicSpecieCatMixture_H

#include "volFields.H"
#include "PtrList.H"
#include "basicCatMixture.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class basicSpecieCatMixture Declaration
\*---------------------------------------------------------------------------*/

class basicSpecieCatMixture
:
    public basicCatMixture
{

protected:

    // Protected data

        //- Table of specie names
        speciesTable species_;

        //- List of specie active flags
        mutable List<bool> active_;

        //- Species mass fractions
        PtrList<volScalarField> Y_;


public:

    //- Run time type information
    TypeName("basicSpecieCatMixture");

    //- The base class of the mixture
    typedef basicSpecieCatMixture basicCatMixtureType;


    // Constructors

        //- Construct from dictionary, species names, mesh and phase name
        basicSpecieCatMixture
        (
            const dictionary&,
            const wordList& specieNames,
            const fvMesh&,
            const word&
        );


    //- Destructor
    virtual ~basicSpecieCatMixture()
    {}


    // Member Functions

        //- Return the table of species
        inline const speciesTable& species() const;

        //- Does the mixture include this specie?
        inline bool contains(const word& specieName) const;

        //- Return true for active species
        inline bool active(label speciei) const;

        //- Return the bool list of active species
        inline const List<bool>& active() const;

        //- Set speciei active
        inline void setActive(label speciei) const;

        //- Set speciei inactive
        inline void setInactive(label speciei) const;

        //- Return the mass-fraction fields
        inline PtrList<volScalarField>& Y();

        //- Return the const mass-fraction fields
        inline const PtrList<volScalarField>& Y() const;

        //- Return the mass-fraction field for a specie given by index
        inline volScalarField& Y(const label i);

        //- Return the const mass-fraction field for a specie given by index
        inline const volScalarField& Y(const label i) const;

        //- Return the mass-fraction field for a specie given by name
        inline volScalarField& Y(const word& specieName);

        //- Return the const mass-fraction field for a specie given by name
        inline const volScalarField& Y(const word& specieName) const;


        // Specie properties

            //- Molecular weight of the given specie [kg/kmol]
            virtual scalar Wi(const label speciei) const = 0;

            //- Enthalpy of formation [J/kg]
            virtual scalar Hf(const label speciei) const = 0;


        // Specie thermo properties

            //- Density [kg/m^3]
            virtual scalar rho
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual scalar Cp
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual scalar HE
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<scalarField> HE
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<volScalarField> HE
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<scalarField> Hs
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField> Hs
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<scalarField> Ha
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual tmp<volScalarField> Ha
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;


        // Specie transport properties

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mu
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Dynamic viscosity [kg/m/s]
            virtual tmp<volScalarField> mu
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Thermal conductivity [W/m/K]
            virtual scalar kappa
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Thermal conductivity [W/m/K]
            virtual tmp<volScalarField> kappa
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Thermal conductivity [W/m/K]
            virtual scalar alphah
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Thermal diffusivity of enthalpy [kg/m/s]
            virtual tmp<volScalarField> alphah
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicSpecieCatMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
