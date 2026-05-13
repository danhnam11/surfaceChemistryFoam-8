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
    Foam::FTMTransport

Description
    Fitting Transport package using Standard Kinetic Theory model. The mixture is 
    calculated based on the modification of the Wilke semi-empirical formulas. 

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Mass Diffusivity is based on mixture averaged model in which binary diffusion 
    coefficients are obtained from Standard Kinetic Theory model.

    It has been validated against NIST and real-fluid based OPPDIF. 
 
    For details of equations in the model, refer to TRANSPORT program 
    in CHEMKIN II @ R. J. Kee et al. CHEMKIN collection, Release 3.6, 
    Reaction Design, Inc., San Diego, CA (2001). 

    Output: 
        Dynamic Viscosity [kg/m.s]
        Thermal Conductivity [W/m.K]    
        Mass Diffusivity [m2/s]    

    developed by:
    Jae Hun Lee - PhD.st
    Advisor: Prof. Chun Sang Yoo 
    Clean combustion & energy research (CCER) Lab - Dept. of Mech. Engineering
    Ulsan Institute of Science and Technology (UNIST) - Korea

SourceFiles
    FTMTransportI.H
    FTMTransport.C

\*---------------------------------------------------------------------------*/

#ifndef FTMTransport_H
#define FTMTransport_H

#include "CoefficientsManager.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class FTMTransport;

template<class Thermo>
inline FTMTransport<Thermo> operator+
(
    const FTMTransport<Thermo>&,
    const FTMTransport<Thermo>&
);

template<class Thermo>
inline FTMTransport<Thermo> operator*
(
    const scalar,
    const FTMTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const FTMTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class FTMTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class FTMTransport
:
    public Thermo
{
    // Private data
    CoefficientsManager internalCoeffs_;
    CoefficientsManager& coeffs_;

    // Private Member Functions

protected:

     //- Members for diffusivity evaluation
	//- list of fitting coefficients

        //- list of mass fraction 
        List<scalar> Ymd_;

        //- list of mole fraction 
        List<scalar> Xmd_;

        List<scalar> wMk_;

      //- protected functions

        //- Temporarily calculate the dynamic viscosity [kg/ms] of pure species
        inline scalar etak
        (
	    label speciei,
            const scalar Wk,
            const scalar p, 
            const scalar T
        ) const;

        //- Temporarily calculate the thermal conductivity [W/mK] of pure species
        inline scalar lamdak
        (
	    label speciei,
            const scalar Wk,
            const scalar p, 
            const scalar T
        ) const;


public:

    // Constructors

        //- Construct from components
        explicit inline FTMTransport
        (
            const Thermo& t,
	    CoefficientsManager& coeffs,
            const List<scalar>& Xmd,
            const List<scalar>& Ymd,
            const List<scalar>& wMk
        );

        //- Construct as named copy
        inline FTMTransport(const word&, const FTMTransport&);

        //- Construct from dictionary
        FTMTransport(const dictionary& dict);

        //- Construct from base thermo and dictionary
       // FTMTransport(const Thermo& t,const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<FTMTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<FTMTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "FTMTransport<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        // inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;

        // 
        // Species diffusivity using Kinetic model for binary diffusion coefficients
        inline scalar Dimix(label speciei, const scalar p, const scalar T) const;

        // Species thermal diffusivity
        inline scalar DimixT(label speciei, const scalar p, const scalar T) const;

        // Binary diffusion coefficient
        inline scalar Dij(label speciei, label speciej, const scalar p, const scalar T) const
        {
            return 1.0;
        }

        // Binary diffusion coefficients Fit function
        inline scalar DijFit
	(
	    label speciei, 
	    label speciej,
	    const scalar p, 
	    const scalar T
	) const;

        //- Update coefficients for mixture transport properties
        inline void updateTRANS
        (
            const List<scalar>& Ymd,
            const List<scalar>& Xmd,
            const List<scalar>& wMk
        );

        //- Update coefficient for Fit transport properties
        inline void updateTRANSFitCoeff
        (
            const List<List<scalar>>& muCoeffsMk,
            const List<List<scalar>>& kappaCoeffsMk,
            const List<List<List<scalar>>>& DijCoeffsMk
        );

        //- Functions return member data for testing mixture  
        // for diffusivity

        inline List<scalar> Ymd() const;
        inline List<scalar> Xmd() const;


        const CoefficientsManager& coeffs() const {return coeffs_;}

        //

    // Member operators

        inline void operator=(const FTMTransport&);

        inline void operator+=(const FTMTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend FTMTransport operator+ <Thermo>
        (
            const FTMTransport&,
            const FTMTransport&
        );

        friend FTMTransport operator* <Thermo>
        (
            const scalar,
            const FTMTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const FTMTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FTMTransportI.H"

#ifdef NoRepository
    #include "FTMTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
