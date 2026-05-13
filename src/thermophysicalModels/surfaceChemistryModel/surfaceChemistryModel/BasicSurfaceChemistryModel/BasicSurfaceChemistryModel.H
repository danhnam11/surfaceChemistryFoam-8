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
    Foam::BasicSurfaceChemistryModel

    developed by Danh Nam Nguyen and Jae Hun Lee,
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).

Description
    Basic chemistry model templated on thermodynamics

SourceFiles
    BasicSurfaceChemistryModelI.H
    BasicSurfaceChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef BasicSurfaceChemistryModel_H
#define BasicSurfaceChemistryModel_H

#include "basicSurfaceChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class BasicSurfaceChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class GasReactionThermo, class SolidReactionThermo>
class BasicSurfaceChemistryModel
:
    public basicSurfaceChemistryModel
{
protected:

    // Protected data

        //- Thermo
        GasReactionThermo& thermo_;
        SolidReactionThermo& catThermo_;

public:

    //- Runtime type information
    TypeName("BasicSurfaceChemistryModel");


    //- Thermo type
    typedef GasReactionThermo reactionThermo;
    typedef SolidReactionThermo reactionCatThermo;

    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        BasicSurfaceChemistryModel,
        thermo,
        (
            GasReactionThermo& thermo,
            SolidReactionThermo& catThermo            
        ),
        (thermo, catThermo)
    );


    // Constructors

        //- Construct from thermo and catThermo
        BasicSurfaceChemistryModel
        (
            GasReactionThermo& thermo,
            SolidReactionThermo& catThermo            
        );

    //- Selector
    static autoPtr
    <
        BasicSurfaceChemistryModel<GasReactionThermo, SolidReactionThermo>
    > New
    (
        GasReactionThermo& thermo,
        SolidReactionThermo& catThermo         
    );

    //- Destructor
    virtual ~BasicSurfaceChemistryModel();


    // Member Functions

        //- Return const access to the thermo package
        inline const GasReactionThermo& thermo() const;

        //- Return non-const access to the thermo package
        inline GasReactionThermo& thermo();

        //- Return const access to the catThermo package
        inline const SolidReactionThermo& catThermo() const;

        //- Return non-const access to the catThermo package
        inline SolidReactionThermo& catThermo();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BasicSurfaceChemistryModel.C"
#endif

#include "BasicSurfaceChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
