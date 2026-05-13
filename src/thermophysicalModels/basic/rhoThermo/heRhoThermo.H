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
    Foam::heRhoThermo

Description
    Energy for a mixture based on density

SourceFiles
    heRhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heRhoThermo_H
#define heRhoThermo_H

#include "rhoThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heRhoThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class heRhoThermo
:
    public heThermo<BasicPsiThermo, MixtureType>
{
    // Private Member Functions

        volScalarField CpNew_; //
        volScalarField kappaNew_; //
        volScalarField WmixNew_; //

        //- Store list of diffusion coefficients
        PtrList<volScalarField> Dimix_;

        //- Store list of thermal diffusion coefficients
        PtrList<volScalarField> DimixT_;

        //- Sotre list of binary diffusion coefficients
        PtrList<PtrList<volScalarField>> Dij_;

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        heRhoThermo(const heRhoThermo<BasicPsiThermo, MixtureType>&);


public:

    //- Runtime type information
    TypeName("heRhoThermo");


    // Constructors

        //- Construct from mesh and phase name
        heRhoThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    //- Destructor
    virtual ~heRhoThermo();


    // Member Functions

        // - Diffusion coefficient of specie ith in the mixture [m^2/s]
        virtual tmp<volScalarField> Dimix( const label speciei) const
        {
            return Dimix_[speciei];
        }

        // - Diffusion coefficient of specie ith in the mixture for patch [m^2/s]
        virtual tmp<scalarField> Dimix( const label speciei, const label patchi) const
        {
            return Dimix_[speciei].boundaryField()[patchi];
        }

        // - Diffusion coefficient of specie ith in the mixture [m^2/s]
        virtual tmp<volScalarField> DimixT( const label speciei) const
        {
            return DimixT_[speciei];
        }

        // - Diffusion coefficient of specie ith in the mixture for patch [m^2/s]
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

        virtual void updateOldField(volScalarField Told, PtrList<volScalarField> Dimix); //Nam

        // Nam
        virtual tmp<volScalarField> CpNew() const
        {
            return CpNew_;
        }
        virtual tmp<volScalarField> kappaNew() const
        {
            return kappaNew_;
        }
        virtual tmp<volScalarField> WmixNew() const
        {
            return WmixNew_;
        }




};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "heRhoThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
