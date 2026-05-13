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
    Foam::PPTransport

Description

    This class is only use for preprocessing of FTM 

    developed by:
    Danh Nam Nguyen and Jae Hun Lee
    Advisor: Prof. Chun Sang Yoo 
    Clean combustion & energy research (CCER) Lab - Dept. of Mech. Engineering
    Ulsan Institute of Science and Technology (UNIST) - Korea

SourceFiles
    PPTransportI.H
    PPTransport.C

\*---------------------------------------------------------------------------*/

#ifndef PPTransport_H
#define PPTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class PPTransport;

template<class Thermo>
inline PPTransport<Thermo> operator+
(
    const PPTransport<Thermo>&,
    const PPTransport<Thermo>&
);

template<class Thermo>
inline PPTransport<Thermo> operator*
(
    const scalar,
    const PPTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const PPTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class PPTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class PPTransport
:
    public Thermo
{
    // Private data

    // Private Member Functions

protected:

     //- Members for diffusivity evaluation
        scalar linearity_;     //linearity
        scalar epsilonOverKb_; //Lennard-Jones well depth
        scalar sigma_;         //Lennard-Jones diameter
        scalar miui_;          //Dipole moment 
        scalar polar_;         //Polarizability
        scalar Zrot_;          //Rotational relaxation number
        scalar Kb_;

        //- list of mass fraction 
        List<scalar> Ymd_;

        //- list of mole fraction 
        List<scalar> Xmd_;

        //- list of epsilonijOverKb[i][j] 
        List<List<scalar>> epsilonijOverKb_;

        //- list of deltaij[i][j] 
        List<List<scalar>> deltaij_;

        //- list of Mij[i][j] 
        List<List<scalar>> Mij_;

        //- list of sigmaij[i][j] 
        List<List<scalar>> sigmaij_;

      //for visc. and cond.
        //- list of all members of k-th specie in the mixture
        List<scalar> linearityMk_;
        List<scalar> epsilonOverKbMk_;
        List<scalar> sigmaMk_;
        List<scalar> miuiMk_;
        List<scalar> polarMk_;
        List<scalar> ZrotMk_;
        List<scalar> wMk_;
        List<List<List<scalar>>> CpCoeffTableMk_;

      //- protected functions
        //- Quadratic interpolation 
        inline scalar qInterpolate(scalar t, List<scalar> X, List<scalar> Y) const;

        // Return collision integral from tablulated data
        //- (Table IV, Monchick and Mason, Journal of Chem. Phys. 35: 1676 (1961))
        inline scalar OMEGA(int N, scalar TR, scalar DR) const;


        //- Funtion return Zrot(T)
        inline scalar ZrotT(scalar Zrot0, scalar epsilonOverKb, scalar T) const;

        //- Temporarily calculate the dynamic viscosity [kg/ms] of pure species
        inline scalar etak
        (
            const scalar Wk,
            const scalar epsilonOverKbk,
            const scalar sigmak,
            const scalar miuik,
            const scalar p, 
            const scalar T
        ) const;

        //- Temporarily calculate the thermal conductivity [W/mK] of pure species
        inline scalar lamdak
        (
            const scalar Wk,
            const scalar linearityk,
            const scalar epsilonOverKbk,
            const scalar sigmak,
            const scalar miuik,
            const scalar Zrotk,
            const List<List<scalar>> CpCoeffTablek,
            const scalar p, 
            const scalar T
        ) const;


public:

    // Constructors

        //- Construct from components
        inline PPTransport
        (
            const Thermo& t,
            const scalar& linearity,
            const scalar& epsilonOverKb,
            const scalar& sigma,
            const scalar& miui,
            const scalar& polar,
            const scalar& Zrot,
            const scalar& Kb,

            const List<scalar>& Ymd,
            const List<scalar>& Xmd,
            const List<List<scalar>>& epsilonijOverKb,
            const List<List<scalar>>& deltaij,
            const List<List<scalar>>& Mij,
            const List<List<scalar>>& sigmaij,

            const List<scalar>& linearityMk,
            const List<scalar>& epsilonOverKbMk,
            const List<scalar>& sigmaMk,
            const List<scalar>& miuiMk,
            const List<scalar>& polarMk,
            const List<scalar>& ZrotMk,
            const List<scalar>& wMk,
            const List<List<List<scalar>>>& CpCoeffTableMk
        );

        //- Construct as named copy
        inline PPTransport(const word&, const PPTransport&);

        //- Construct from dictionary
        PPTransport(const dictionary& dict);

        //- Construct from base thermo and dictionary
       // PPTransport(const Thermo& t,const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<PPTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<PPTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "preProcessingFTMTransport<" + Thermo::typeName() + '>';
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

        // Species thermal diffusivity using Kinetic model for binary diffusion coefficients
        inline scalar DimixT(label speciei, const scalar p, const scalar T) const;


        // Binary diffusion coefficient
        inline scalar Dij
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
            const List<List<scalar>>& epsilonijOverKb,
            const List<List<scalar>>& deltaij,
            const List<List<scalar>>& Mij,
            const List<List<scalar>>& sigmaij,
            
            const List<scalar>& linearityMk,
            const List<scalar>& epsilonOverKbMk,
            const List<scalar>& sigmaMk,
            const List<scalar>& miuiMk,
            const List<scalar>& polarMk,
            const List<scalar>& ZrotMk,
            const List<scalar>& wMk,
            const List<List<List<scalar>>>& CpCoeffTableMk 
        );

        //- Functions return member data for testing mixture  
        // for diffusivity
        inline scalar linearity() const;
        inline scalar epsilonOverKb() const;
        inline scalar sigma() const;
        inline scalar miui() const;
        inline scalar polar() const;
        inline scalar Zrot() const;
        inline scalar Kb() const;

        inline List<scalar> Ymd() const;
        inline List<scalar> Xmd() const;
        inline List<List<scalar>> epsilonijOverKb() const;
        inline List<List<scalar>> deltaij() const;
        inline List<List<scalar>> Mij() const;
        inline List<List<scalar>> sigmaij() const;

        //

    // Member operators

        inline void operator=(const PPTransport&);

        inline void operator+=(const PPTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend PPTransport operator+ <Thermo>
        (
            const PPTransport&,
            const PPTransport&
        );

        friend PPTransport operator* <Thermo>
        (
            const scalar,
            const PPTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const PPTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PPTransportI.H"

#ifdef NoRepository
    #include "PPTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
