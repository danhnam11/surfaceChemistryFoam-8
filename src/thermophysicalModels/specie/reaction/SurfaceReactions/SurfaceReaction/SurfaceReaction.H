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
    Foam::SurfaceReaction

    developed by Danh Nam Nguyen,
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).

Description
    Simple extension of SolidThermoType to handle reaction kinetics in addition
    to the equilibrium thermodynamics already handled.

SourceFiles
    SurfaceReactionI.H
    SurfaceReaction.C

\*---------------------------------------------------------------------------*/

#ifndef SurfaceReaction_H
#define SurfaceReaction_H

#include "speciesTable.H"
#include "catSpecieCoeffs.H"
#include "HashPtrTable.H"
#include "scalarField.H"
#include "simpleMatrix.H"
#include "Tuple2.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class GasThermoType, class SolidThermoType>
class SurfaceReaction;

template<class GasThermoType, class SolidThermoType>
inline Ostream& operator<<(Ostream&, const SurfaceReaction<GasThermoType, SolidThermoType>&);

class objectRegistry;

/*---------------------------------------------------------------------------*\
                           Class SurfaceReaction Declaration
\*---------------------------------------------------------------------------*/

template<class GasThermoType, class SolidThermoType>
class SurfaceReaction
:
    public SolidThermoType::thermoType
{

public:

    // Static data

        //- Number of un-named reactions
        static label nUnNamedSurfaceReactions;

        //- Default temperature limits of applicability of reaction rates
        static scalar TlowDefault, ThighDefault;

private:

    // Private Data

        //- Name of reaction
        const word name_;

        //- List of specie names present in reaction system
        //const speciesTable& species_;

        //- Temperature limits of applicability of reaction rates
        scalar Tlow_, Thigh_;

        //- Specie info for the left-hand-side of the reaction
        List<catSpecieCoeffs> lhs_;

        //- Specie info for the right-hand-side of the reaction
        List<catSpecieCoeffs> rhs_;

        //- gas specie info for the left-hand-side of the reaction
        List<catSpecieCoeffs> lhsG_;

        //- gas specie info for the right-hand-side of the reaction
        List<catSpecieCoeffs> rhsG_;

    // Private Member Functions

        //- Construct reaction thermo
        void setThermo(const HashPtrTable<SolidThermoType>& thermoDatabase);

        //- Return new reaction ID for un-named reactions
        label getNewSurfaceReactionID();

protected:
        //- List of specie names present in reaction system
        const speciesTable& gasSpecies_;    
        const speciesTable& species_;

        label nG_;
        label nGasSpecies_;
        label nSolidSpecies_;
        List<scalar> gasWk_;


public:

    //- Runtime type information
    TypeName("SurfaceReaction");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            SurfaceReaction,
            dictionary,
            (
                const speciesTable& gasSpecies,
                const speciesTable& species,                
                const HashPtrTable<GasThermoType>& gasThermoDatabase,
                const HashPtrTable<SolidThermoType>& thermoDatabase,
                const dictionary& dict
            ),
            (gasSpecies, species, gasThermoDatabase, thermoDatabase, dict)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            SurfaceReaction,
            objectRegistry,
            (
                const speciesTable& gasSpecies,
                const speciesTable& species,                
                const HashPtrTable<GasThermoType>& gasThermoDatabase,
                const HashPtrTable<SolidThermoType>& thermoDatabase,                
                const objectRegistry& ob,
                const dictionary& dict
            ),
            (gasSpecies, species, gasThermoDatabase, thermoDatabase, ob, dict)
        );


    // Constructors

        /*
        //- Construct from components
        SurfaceReaction
        (
            const speciesTable& gasSpecies,
            const speciesTable& species,
            const List<catSpecieCoeffs>& lhsG,
            const List<catSpecieCoeffs>& rhsG,                        
            const List<catSpecieCoeffs>& lhs,
            const List<catSpecieCoeffs>& rhs,
            const HashPtrTable<GasThermoType>& gasThermoDatabase            
            const HashPtrTable<SolidThermoType>& thermoDatabase
        );
        */

        //- Construct as copy given new speciesTable
        SurfaceReaction
        (
            const SurfaceReaction<GasThermoType, SolidThermoType>&, 
            const speciesTable& gasSpecies,
            const speciesTable& species            
        );


        //- Construct from dictionary
        SurfaceReaction
        (
            const speciesTable& gasSpecies,
            const speciesTable& species,                
            const HashPtrTable<GasThermoType>& gasThermoDatabase,
            const HashPtrTable<SolidThermoType>& thermoDatabase,  
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<SurfaceReaction<GasThermoType, SolidThermoType>> clone() const = 0;

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<SurfaceReaction<GasThermoType, SolidThermoType>> clone
        (
            const speciesTable& gasSpecies,
            const speciesTable& species            
        ) const = 0;


    // Selectors

        //- Return a pointer to new patchField created on freestore
        //  from dictionary
        static autoPtr<SurfaceReaction<GasThermoType, SolidThermoType>> New
        (
            const speciesTable& gasSpecies,
            const speciesTable& species,                
            const HashPtrTable<GasThermoType>& gasThermoDatabase,
            const HashPtrTable<SolidThermoType>& thermoDatabase,  
            const dictionary& dict
        );

        //- Return a pointer to new patchField created on freestore
        //  from objectRegistry and dictionary
        static autoPtr<SurfaceReaction<GasThermoType, SolidThermoType>> New
        (
            const speciesTable& gasSpecies,
            const speciesTable& species,                
            const HashPtrTable<GasThermoType>& gasThermoDatabase,
            const HashPtrTable<SolidThermoType>& thermoDatabase,  
            const objectRegistry& ob,
            const dictionary& dict
        );

        //- Return a pointer to new patchField created on freestore from dict.
        //  This is used to construct a single reaction given a speciesThermo
        //  list from which a temporary thermoDatabase is created.
        static autoPtr<SurfaceReaction<GasThermoType, SolidThermoType>> New
        (
            const speciesTable& gasSpecies,
            const speciesTable& species,
            const PtrList<GasThermoType>& gasSpeciesThermo,
            const PtrList<SolidThermoType>& speciesThermo,            
            const dictionary& dict
        );


    //- Destructor
    virtual ~SurfaceReaction()
    {}


    // Member Functions

        // Access

            //- Return the name of the reaction
            inline const word& name() const;

            //- Return the lower temperature limit for the reaction
            inline scalar Tlow() const;

            //- Return the upper temperature limit for the reaction
            inline scalar Thigh() const;

            //- Return the components of the left hand side
            inline const List<catSpecieCoeffs>& lhs() const;

            //- Return the components of the right hand side
            inline const List<catSpecieCoeffs>& rhs() const;

            //- Return the gas components of the left hand side
            inline const List<catSpecieCoeffs>& lhsG() const;

            //- Return the gas components of the right hand side
            inline const List<catSpecieCoeffs>& rhsG() const;

            //- Return the specie list
            const speciesTable& species() const;

            //- Return the gas specie list
            const speciesTable& gasSpecies() const;


        // SurfaceReaction rate coefficients

            //- Forward reaction rate
            void ddot
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarField& d
            ) const;

            //- Backward reaction rate
            void fdot
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarField& f
            ) const;

            //- Net reaction rate for individual species
            void omega
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarField& dcdt,
                List<scalar> siteDensityList,
                List<label> siteNumberList                
            ) const;

            //- Net reaction rate
            scalar omega
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalar& pf,
                scalar& cf,
                label& lRef,
                scalar& pr,
                scalar& cr,
                label& rRef,
                List<scalar> siteDensityList,
                List<label> siteNumberList                
            ) const;

        // SurfaceReaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                List<scalar> siteDensityList,
                List<label> siteNumberList                
            ) const = 0;

            //- Reverse rate constant from the given forward rate constant
            virtual scalar kr
            (
                const scalar kfwd,
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const = 0;

            //- Reverse rate constant
            virtual scalar kr
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const = 0;


        // Jacobian coefficients

            //- Derivative of the net reaction rate for each species involved
            //  w.r.t. the species concentration
            void dwdc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarSquareMatrix& J,
                scalarField& dcdt,
                scalar& omegaI,
                scalar& kfwd,
                scalar& kbwd,
                const bool reduced,
                const List<label>& c2s,
                List<scalar> siteDensityList,
                List<label> siteNumberList                
            ) const;

            //- Derivative of the net reaction rate for each species involved
            //  w.r.t. the temperature
            void dwdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                const scalar omegaI,
                const scalar kfwd,
                const scalar kbwd,
                scalarSquareMatrix& J,
                const bool reduced,
                const List<label>& c2s,
                const label indexT
            ) const;

            //- Temperature derivative of forward rate
            virtual scalar dkfdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const = 0;

            //- Temperature derivative of reverse rate
            virtual scalar dkrdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                const scalar dkfdT,
                const scalar kr
            ) const = 0;

            //- Third-body efficiencies (beta = 1-alpha)
            //  non-empty only for third-body reactions
            //  with enhanced molecularity (alpha != 1)
            // virtual const List<Tuple2<label, scalar>>& beta() const = 0;

           //- Coverage modification coeffs
           //  non-empty only for surface coverage dependent reactions
            virtual const List<Tuple2<label, scalarField>>& coverageModificationCoeffs() const = 0;


            //- Species concentration derivative of the pressure dependent term
            virtual void dcidc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarField& dcidc
            ) const = 0;

            //- Temperature derivative of the pressure dependent term
            virtual scalar dcidT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const = 0;


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SurfaceReaction<GasThermoType, SolidThermoType>&) = delete;


    // Ostream Operator

        friend Ostream& operator<< <GasThermoType, SolidThermoType>
        (
            Ostream&,
            const SurfaceReaction<GasThermoType, SolidThermoType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SurfaceReactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SurfaceReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
