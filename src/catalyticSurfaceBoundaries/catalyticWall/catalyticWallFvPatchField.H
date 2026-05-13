/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::catalyticWallFvPatchField

    developed by Danh Nam Nguyen
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).

Description
    This boundary condition is only used for site species. 
    The value of site species fraction has been updated directly inside 
    the surface chemistry library.

Usage

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type         catalyticWall;
        value        uniform 0;  // value is read as initial value
                                 // before solving surface reactions
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField
    Foam::Function1s

SourceFiles
    catalyticWallFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef catalyticWallFvPatchField_H
#define catalyticWallFvPatchField_H

#include "fixedValueFvPatchFields.H"
//#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class catalyticWallFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class catalyticWallFvPatchField
:
    public fixedValueFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("catalyticWall");


    // Constructors

        //- Construct from patch and internal field
        catalyticWallFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        catalyticWallFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given catalyticWallFvPatchField
        //  onto a new patch
        catalyticWallFvPatchField
        (
            const catalyticWallFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        catalyticWallFvPatchField
        (
            const catalyticWallFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new catalyticWallFvPatchField<Type>(*this)
            );
        }

        //- Copy constructor setting internal field reference
        catalyticWallFvPatchField
        (
            const catalyticWallFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new catalyticWallFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "catalyticWallFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
