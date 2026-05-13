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
    Foam::combustionModels::noCombustion

Description
    Dummy combustion model for 'no combustion'

SourceFiles
    noCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef noCombustion_H
#define noCombustion_H

#include "ThermoCombustion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                        Class noCombustion Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class noCombustion
:
    public ThermoCombustion<ReactionThermo>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noCombustion
        (
            const word& modelType,
            const ReactionThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        noCombustion(const noCombustion&);


    //- Destructor
    virtual ~noCombustion();


    // Member Functions

        //- Correct combustion rate
        virtual void correct();

        //- Fuel consumption rate matrix
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s^3]
        virtual tmp<volScalarField> Qdot() const;

        //- Update properties from given dictionary
        virtual bool read();

        //- Update Mass Fractions // JH for splitting algorithm
        virtual void updateMassFractions()
	{
	    Info << "updateMassFractions at the noCombustion " << endl;
	};

        //virtual void updateTfromTchem(volScalarField& T)
        virtual void updateYiHe(PtrList<volScalarField>& Y, volScalarField& he)
	{
	    Info << " Do nothing here at updateTfromTchem(T) func. " << endl;
	};

        virtual void solveChemistry
        (
            PtrList<volScalarField> Yold,
            volScalarField p,
            volScalarField T,
            PtrList<volScalarField>& Ynew,
            volScalarField& Qdot
        )
        {
            Info << " Do nothing here at solveChemistry func. " << endl;
        };


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const noCombustion&) = delete;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "noCombustion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
