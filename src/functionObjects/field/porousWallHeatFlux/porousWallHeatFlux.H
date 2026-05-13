/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::functionObjects::porousWallHeatFlux

Description
    Calculates and write the heat-flux at wall patches as the
    volScalarField field 'porousWallHeatFlux'.

    All wall patches are included by default; to restrict the calculation to
    certain patches, use the optional 'patches' entry.

    Example of function object specification:
    \verbatim
    porousWallHeatFlux1
    {
        type        porousWallHeatFlux;
        libs        ("libfieldFunctionObjects.so");
        ...
        region      fluid;
        patches     (".*Wall");
    }
    \endverbatim

Usage
    \table
        Property | Description                | Required   | Default value
        type     | type name: porousWallHeatFlux    | yes        |
        patches  | list of patches to process | no         | all wall patches
        region   | region to be evaluated     | no         | default region
    \endtable

Note
    Writing field 'porousWallHeatFlux' is done by default, but it can be overridden by
    defining an empty \c objects list. For details see writeLocalObjects.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    porousWallHeatFlux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_porousWallHeatFlux_H
#define functionObjects_porousWallHeatFlux_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
#include "HashSet.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class porousWallHeatFlux Declaration
\*---------------------------------------------------------------------------*/

class porousWallHeatFlux
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{

protected:

    // Protected data

        //- Optional list of patches to process
        labelHashSet patchSet_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate the heat-flux
        tmp<volScalarField> calcWallHeatFlux
        (
            const volScalarField& alpha,
            const volScalarField& he,
            const volScalarField& kappaEff,
            const volScalarField& Cp
        );


public:

    //- Runtime type information
    TypeName("porousWallHeatFlux");


    // Constructors

        //- Construct from Time and dictionary
        porousWallHeatFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        porousWallHeatFlux(const porousWallHeatFlux&) = delete;


    //- Destructor
    virtual ~porousWallHeatFlux();


    // Member Functions

        //- Read the porousWallHeatFlux data
        virtual bool read(const dictionary&);

        //- Calculate the wall heat-flux
        virtual bool execute();

        //- Write the wall heat-flux
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const porousWallHeatFlux&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
