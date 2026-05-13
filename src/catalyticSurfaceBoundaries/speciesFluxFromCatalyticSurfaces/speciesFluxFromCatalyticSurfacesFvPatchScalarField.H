/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::speciesFluxFromCatalyticSurfacesFvPatchScalarField

    developed by Danh Nam Nguyen and Jae Hun Lee,
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).

Description

    Species flux model for catalytic surfaces (boundaries). 
    Effect of stefan flow has been added in the surface chemistry library.

Usage

    \verbatim
    boundaryField
    {
        <patchName>    // e.g., catalyticWalls, change if the patch name is different
        {
            type        speciesFluxFromCatalyticSurfaces;
            gradient    uniform 0.0;
            value       $internalField;
        }
        ...
    }
    \endverbatim

See also
    Foam::fixedGradientFvPatchField

SourceFiles
    speciesFluxFromCatalyticSurfacesFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef speciesFluxFromCatalyticSurfacesFvPatchScalarField_H
#define speciesFluxFromCatalyticSurfacesFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class speciesFluxFromCatalyticSurfacesFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class speciesFluxFromCatalyticSurfacesFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{


public:

    //- Runtime type information
    TypeName("speciesFluxFromCatalyticSurfaces");


    // Constructors

        //- Construct from patch and internal field
        speciesFluxFromCatalyticSurfacesFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        speciesFluxFromCatalyticSurfacesFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given speciesFluxFromCatalyticSurfacesFvPatchScalarField
        // onto a new patch
        speciesFluxFromCatalyticSurfacesFvPatchScalarField
        (
            const speciesFluxFromCatalyticSurfacesFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        speciesFluxFromCatalyticSurfacesFvPatchScalarField
        (
            const speciesFluxFromCatalyticSurfacesFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new speciesFluxFromCatalyticSurfacesFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        speciesFluxFromCatalyticSurfacesFvPatchScalarField
        (
            const speciesFluxFromCatalyticSurfacesFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new speciesFluxFromCatalyticSurfacesFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
