/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::stefanFlowVelocityFvPatchVectorField

    developed by Jae Hun Lee,
    Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering,
    Ulsan National Institute of Science and Technology (UNIST), Korea 
    (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/).

Description

    stefan flow velocity boundary condition based on 
    Surface-normal fixed value vector boundary condition 

Usage

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            stefanFlowVelocity;
        refValue        uniform 0; // dummy value, doest affect to BCs
    }
    \endverbatim

Note
    Sign conventions:
    - the value is positive for outward-pointing vectors

See also
    Foam::fixedValueFvPatchField

SourceFiles
    stefanFlowVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef stefanFlowVelocityFvPatchVectorField_H
#define stefanFlowVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class stefanFlowVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class stefanFlowVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        scalarField refValue_;


public:

    //- Runtime type information
    TypeName("stefanFlowVelocity");


    // Constructors

        //- Construct from patch and internal field
        stefanFlowVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        stefanFlowVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  stefanFlowVelocityFvPatchVectorField
        //  onto a new patch
        stefanFlowVelocityFvPatchVectorField
        (
            const stefanFlowVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        stefanFlowVelocityFvPatchVectorField
        (
            const stefanFlowVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new stefanFlowVelocityFvPatchVectorField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        stefanFlowVelocityFvPatchVectorField
        (
            const stefanFlowVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new stefanFlowVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchVectorField&, const labelList&);


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
