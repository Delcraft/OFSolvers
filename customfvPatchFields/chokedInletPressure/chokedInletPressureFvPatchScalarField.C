/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "chokedInletPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chokedInletPressureFvPatchScalarField::chokedInletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_("phi"),
    p0_(1e5),
    T0_(300.0),
    gamma_(1.4),
    R_(287.04),
    Cd_(1.0)
{
    // Initialize mixedFvPathField members
    this->refValue() = *this;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}

Foam::chokedInletPressureFvPatchScalarField::chokedInletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    p0_(readScalar(dict.lookup("p0"))),
    T0_(readScalar(dict.lookup("T0"))),
    gamma_(readScalar(dict.lookup("gamma"))),
    R_(readScalar(dict.lookup("R"))),
    Cd_(readScalar(dict.lookup("Cd")))
{
    this->refValue() = Field<scalar>("inletValue", dict, p.size());

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            Field<scalar>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(this->refValue());
    }

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}

Foam::chokedInletPressureFvPatchScalarField::chokedInletPressureFvPatchScalarField
(
    const chokedInletPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    p0_(ptf.p0_),
    T0_(ptf.T0_),
    gamma_(ptf.gamma_),
    R_(ptf.R_),
    Cd_(ptf.Cd_)
{}

Foam::chokedInletPressureFvPatchScalarField::chokedInletPressureFvPatchScalarField
(
    const chokedInletPressureFvPatchScalarField& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    phiName_(ptf.phiName_),
    p0_(ptf.p0_),
    T0_(ptf.T0_),
    gamma_(ptf.gamma_),
    R_(ptf.R_),
    Cd_(ptf.Cd_)
{}

Foam::chokedInletPressureFvPatchScalarField::chokedInletPressureFvPatchScalarField
(
    const chokedInletPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    phiName_(ptf.phiName_),
    p0_(ptf.p0_),
    T0_(ptf.T0_),
    gamma_(ptf.gamma_),
    R_(ptf.R_),
    Cd_(ptf.Cd_)
{}

// ADD CONSTRUCTORS HERE

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::chokedInletPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    Info << "Came to updateCoeffs()" << endl;

    // Patch properties

    // Internal pressure field
    
    const fvPatchField<scalar>& p =
        db().lookupObject<volScalarField>
        (
            internalField().name()
        ).boundaryField()[patch().index()];
    

    const fvPatchField<scalar>& rho_ = 
        db().lookupObject<volScalarField>
        (
            "rho"
        ).boundaryField()[patch().index()];
/*
    const fvMesh& mesh =
        db().lookupObject<volScalarField>
        (
            internalField().name()
        ).mesh();
*/
    scalar gammaM1 = gamma_ - 1.0;

    // Compute mean pressure
    scalar pmean = gSum(p*this->patch().magSf())/gSum(this->patch().magSf());

    scalar rhomean = gSum(rho_*this->patch().magSf())/gSum(this->patch().magSf());
/*
    scalar totalVolume = scalar(0.0);

    scalar pmeanCell = 0.0;

    const scalarField& pInt_ = db().lookupObject<volScalarField>("p").internalField();

    forAll(patch().faceCells(),cellI)
    {
//        Info << "pCell: " << pInt_[cellI] << endl;
        pmeanCell += pInt_[cellI]*mesh.V()[cellI];
        totalVolume += mesh.V()[cellI];
    }

    pmeanCell = pmeanCell/totalVolume;
*/

    Info << endl;
    Info << "pmean: " << pmean << endl;
    
    // Compute injection velocity

    scalar massflow = 0.0;
    scalar u_inj = 0.0;

    if (pmean/p0_ > 1)
    {
        FatalErrorInFunction
            << "pmean/p0_ > 1 : backflow occured!"
            << "\n    on patch " << patch().name()
            << exit(FatalError);
    }
    else
    {
        massflow = Cd_*p0_/sqrt(R_*T0_)*pow(pmean/p0_,1.0/gamma_)
            *sqrt(2.0*gamma_/gammaM1
            *(1.0-pow(pmean/p0_,gammaM1/gamma_)));

        u_inj = massflow/rhomean;
    }

    Info << "u_inj: " << u_inj << endl;

    /*
    const fvPatchField<scalar>& rho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);
    */

    // Injection Mach number
    scalar Mach = u_inj / sqrt(gamma_*R_*T0_);
    Info << "Mach: " << Mach << endl;

    if (Mach < 1)
    {
        // Subsonic regime, pressure is taken from internal field
        // w = 0.0, apply zero gradient condition
        Info << "Subsonic regime" << endl;
        this->valueFraction() = 0.0;
        this->refGrad() = 0.0;
    }
    else
    {
        // Sonic regime, pressure is imposed at the boundary
        // w = 0, refValue = static pressure at injection surface
        Info << "Sonic regime" << endl;
        this->refValue() = pow(2.0/(gamma_ + 1),gamma_/gammaM1)*p0_;
        this->valueFraction() = 1.0;
    }

    mixedFvPatchField<scalar>::updateCoeffs();
}

void Foam::chokedInletPressureFvPatchScalarField::write(Ostream& os) const
{

    Info << "write function" << endl;
    fvPatchField<scalar>::write(os);
    os.writeKeyword("p0") << p0_
        << token::END_STATEMENT << nl;
    os.writeKeyword("T0") << T0_
        << token::END_STATEMENT << nl;
    os.writeKeyword("gamma") << gamma_
        << token::END_STATEMENT << nl;
    os.writeKeyword("R") << R_
        << token::END_STATEMENT << nl;
    os.writeKeyword("Cd") << Cd_
        << token::END_STATEMENT << nl;
    this->refValue().writeEntry("inletValue", os);
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        chokedInletPressureFvPatchScalarField
    );
}

// ************************************************************************* //