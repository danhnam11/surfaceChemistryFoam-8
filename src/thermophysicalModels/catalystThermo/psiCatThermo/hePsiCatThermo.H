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

Class
    Foam::hePsiCatThermo

    developed by Danh Nam Nguyen
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).


Description
    Energy for a mixture based on compressibility

SourceFiles
    hePsiCatThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hePsiCatThermo_H
#define hePsiCatThermo_H

#include "psiCatThermo.H"
#include "heCatThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class hePsiCatThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class hePsiCatThermo
:
    public heCatThermo<BasicPsiThermo, MixtureType>
{

	// const fvMesh& mesh_;

    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        //Nam
        label numberOfSpecies_;

        //- Store list of diffusion coefficients
        PtrList<volScalarField> Dimix_;

        //- Store list of thermal diffusion coefficients
        PtrList<volScalarField> DimixT_;

        //- Store list of binary diffusion coefficients
        List<PtrList<volScalarField>> Dij_;

        const PtrList<typename MixtureType::thermoType>& specieThermos_;

	    const fvMesh& mesh_;

        //- Construct as copy (not implemented)
        hePsiCatThermo(const hePsiCatThermo<BasicPsiThermo, MixtureType>&);

public:

    //- Runtime type information
    TypeName("hePsiCatThermo");


    // Constructors

        //- Construct from mesh and phase name
        hePsiCatThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    //- Destructor
    virtual ~hePsiCatThermo();


    // Member Functions

        // - Diffusion coefficient of specie ith in the mixture [m^2/s] - Nam
        virtual tmp<volScalarField> Dimix( const label speciei) const
        {
            return Dimix_[speciei];
        }

        // - Diffusion coefficient of specie ith in the mixture for patch [m^2/s] - Nam
        virtual tmp<scalarField> Dimix( const label speciei, const label patchi) const
        {
            return Dimix_[speciei].boundaryField()[patchi];
        }

        // - Thermal diffusion coefficient of specie ith in the mixture [m^2/s] - Nam
        virtual tmp<volScalarField> DimixT( const label speciei) const
        {
            return DimixT_[speciei];
        }

        // - Thermal diffusion coefficient of specie ith in the mixture for patch [m^2/s] - Nam
        virtual tmp<scalarField> DimixT( const label speciei, const label patchi) const
        {
            return DimixT_[speciei].boundaryField()[patchi];
        }

        // - Binary diffusion coefficient of species ith and species jth in the mixture [m^2/s] - JH
        virtual tmp<volScalarField> Dij( const label speciei, const label speciej) const
        {
            return Dij_[speciei][speciej];
        }

        // - Binary diffusion coefficient of species ith and species jth in the mixture for patch [m^2/s] - JH
        virtual tmp<scalarField> Dij( const label speciei, const label speciej,  const label patchi) const
        {
            return Dij_[speciei][speciej].boundaryField()[patchi];
        }

        //- Update properties
        virtual void correct();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "hePsiCatThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
