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
    Foam::hePsiThermo

Description
    Energy for a mixture based on compressibility

SourceFiles
    hePsiThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hePsiThermo_H
#define hePsiThermo_H

#include "psiThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class hePsiThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class hePsiThermo
:
    public heThermo<BasicPsiThermo, MixtureType>
{

	// const fvMesh& mesh_;

    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        //- Initialize the thermo variables, Nam
        void initialize();

        //- Calculate the thermo variables 
        void calculateUsingCoTHERMOnlyT();

        //- Calculate the thermo variables
        void calculateUsingCoTHERM();

        //- Calculate for the transport properties pre process -JH
        void calculateTransportPreProcessing();

        //Nam
        label numberOfSpecies_;

        //- Store list of diffusion coefficients
        PtrList<volScalarField> Dimix_;

        //- Store list of thermal diffusion coefficients
        PtrList<volScalarField> DimixT_;

        //- Store list of binary diffusion coefficients
        List<PtrList<volScalarField>> Dij_;

        volScalarField CpNew_;       
        volScalarField kappaNew_;       
        volScalarField WmixNew_;       
        const PtrList<typename MixtureType::thermoType>& specieThermos_;

	    const fvMesh& mesh_;

        dictionary thermoDict_;
        bool usingDetailedTransportModel_;
        //bool usingTransportPreProcessing_;
        bool usingPreProcessingFTM_;
        bool usingCoTHERM_;
        bool usingCoTHERMOnlyT_;        
        scalar epsilonOnlyT_;                
        scalar epsilonT_;
        scalar epsilonP_;        
        volScalarField residualT_;
        volScalarField residualP_;            
        volScalarField flagCell_; 

        // to avoid coTHERM copy forever 
        // Number of consecutive CoTHERM-based reuse steps per cell
        volScalarField coTHERMStepCount_;
        // Maximum allowed reuse count before forcing detailed recalculation
        scalar maxCoTHERMStepCount_;

        //- coTHERM modes
        int mode_;

        //- Construct as copy (not implemented)
        hePsiThermo(const hePsiThermo<BasicPsiThermo, MixtureType>&);

public:

    //- Runtime type information
    TypeName("hePsiThermo");


    // Constructors

        //- Construct from mesh and phase name
        hePsiThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    //- Destructor
    virtual ~hePsiThermo();


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
        virtual void updateOldField(volScalarField Told, PtrList<volScalarField> Dimix); //Nam

        // - Nam
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
    #include "hePsiThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
