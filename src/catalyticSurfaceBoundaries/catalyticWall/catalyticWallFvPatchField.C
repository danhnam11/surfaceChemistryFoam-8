/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "catalyticWallFvPatchField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::catalyticWallFvPatchField<Type>::catalyticWallFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::catalyticWallFvPatchField<Type>::catalyticWallFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::catalyticWallFvPatchField<Type>::catalyticWallFvPatchField
(
    const catalyticWallFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::catalyticWallFvPatchField<Type>::catalyticWallFvPatchField
(
    const catalyticWallFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::catalyticWallFvPatchField<Type>::catalyticWallFvPatchField
(
    const catalyticWallFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::catalyticWallFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //const scalar t = this->db().time().timeOutputValue();
    //Type meanValue = meanValue_->value(t);

    //Field<Type> newValues(this->patchInternalField());
    // Field<Type> newValues(*this);    

    /*
    Type meanValuePsi =
        gSum(this->patch().magSf()*newValues)
       /gSum(this->patch().magSf());

    if (mag(meanValue) > small && mag(meanValuePsi)/mag(meanValue) > 0.5)
    {
        newValues *= mag(meanValue)/mag(meanValuePsi);
    }
    else
    {
        newValues += (meanValue - meanValuePsi);
    }
    */

    // this->operator==(newValues);

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::catalyticWallFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
