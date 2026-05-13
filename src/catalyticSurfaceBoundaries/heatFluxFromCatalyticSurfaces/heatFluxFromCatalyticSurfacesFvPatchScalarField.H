/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::heatFluxFromCatalyticSurfacesFvPatchScalarField

    developed by Danh Nam Nguyen and Jae Hun Lee,
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).

Description

    Heat flux model for catalytic surfaces (boundaries)

Usage
    \verbatim
    boundaryField
    {
        <patchName>    // e.g., catalyticWalls, change if the patch name is different
        {
            type        heatFluxFromCatalyticSurfaces;
            gradient    uniform 0.0;
            value       $internalField;
        }
        ...
    }
    \endverbatim

See also
    Foam::fixedGradientFvPatchField

SourceFiles
    heatFluxFromCatalyticSurfacesFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef heatFluxFromCatalyticSurfacesFvPatchScalarField_H
#define heatFluxFromCatalyticSurfacesFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class heatFluxFromCatalyticSurfacesFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class heatFluxFromCatalyticSurfacesFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("heatFluxFromCatalyticSurfaces");


    // Constructors

        //- Construct from patch and internal field
        heatFluxFromCatalyticSurfacesFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        heatFluxFromCatalyticSurfacesFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given heatFluxFromCatalyticSurfacesFvPatchScalarField
        // onto a new patch
        heatFluxFromCatalyticSurfacesFvPatchScalarField
        (
            const heatFluxFromCatalyticSurfacesFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        heatFluxFromCatalyticSurfacesFvPatchScalarField
        (
            const heatFluxFromCatalyticSurfacesFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new heatFluxFromCatalyticSurfacesFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        heatFluxFromCatalyticSurfacesFvPatchScalarField
        (
            const heatFluxFromCatalyticSurfacesFvPatchScalarField&,
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
                new heatFluxFromCatalyticSurfacesFvPatchScalarField(*this, iF)
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
