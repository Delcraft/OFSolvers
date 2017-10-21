/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::chokedInletPressureFvPatchScalarField

Group
    grpInletBoundaryConditions

SourceFiles
    chokedInletPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef chokedInletPressureFvPatchScalarField_H
#define chokedInletPressureFvPatchScalarField_H

#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class chokedInletPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class chokedInletPressureFvPatchScalarField
:
    public mixedFvPatchField<scalar>
{
    // Private data

    word phiName_;

    // Total pressure
    scalar p0_;

    // Total temperature
    scalar T0_;

    // Ratio of specific heats
    scalar gamma_;

    // Specific gas constant
    scalar R_;

    // Inlet discharge coefficient
    scalar Cd_;

public:

    //- Runtime type information
    TypeName("chokedInletPressure");


    // Constructors
        //- Construct from patch and internal field
        chokedInletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        chokedInletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given chokedInletPressureFvPatchScalarField
        //  onto a new patch
        chokedInletPressureFvPatchScalarField
        (
            const chokedInletPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        chokedInletPressureFvPatchScalarField
        (
            const chokedInletPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new chokedInletPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        chokedInletPressureFvPatchScalarField
        (
            const chokedInletPressureFvPatchScalarField&,
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
                new chokedInletPressureFvPatchScalarField(*this, iF)
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