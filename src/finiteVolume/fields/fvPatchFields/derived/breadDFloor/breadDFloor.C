/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "breadDFloor.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// template<class Type>
bool Foam::breadDFloorFvPatchVectorField::readMixedEntries
(
    const dictionary& dict,
    IOobjectOption::readOption readOpt
)
{
    if (!IOobjectOption::isAnyRead(readOpt)) return false;
    const auto& p = fvPatchFieldBase::patch();


    // If there is a 'refValue', also require all others
    const auto* hasValue = dict.findEntry("refValue", keyType::LITERAL);

    if (!hasValue && IOobjectOption::isReadOptional(readOpt))
    {
        return false;
    }

    const auto* hasGrad = dict.findEntry("refGradient", keyType::LITERAL);
    const auto* hasFrac = dict.findEntry("valueFraction", keyType::LITERAL);

    // Combined error message on failure
    if (!hasValue || !hasGrad || !hasFrac)
    {
        FatalIOErrorInFunction(dict)
            << "Required entries:";

        if (!hasValue) FatalIOError << " 'refValue'";
        if (!hasGrad)  FatalIOError << " 'refGradient'";
        if (!hasFrac)  FatalIOError << " 'valueFraction'";

        FatalIOError
            << " : missing for patch " << p.name()
            << " : in dictionary " << dict.relativeName() << nl
            << exit(FatalIOError);
    }

    // Everything verified - can assign
    refValue_.assign(*hasValue, p.size());
    refGrad_.assign(*hasGrad, p.size());
    valueFraction_.assign(*hasFrac, p.size());
    dict.readEntry("floorPos", floorPos_);
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// template<class Type>
Foam::breadDFloorFvPatchVectorField::breadDFloorFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fvPatchVectorField(p, iF),
    refValue_(p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size()),
    source_(p.size(), Zero)
{}


// template<class Type>
Foam::breadDFloorFvPatchVectorField::breadDFloorFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const Foam::zero
)
:
    fvPatchVectorField(p, iF),
    refValue_(p.size(), Zero),
    refGrad_(p.size(), Zero),
    valueFraction_(p.size(), Zero),
    source_(p.size(), Zero)
{}


// template<class Type>
Foam::breadDFloorFvPatchVectorField::breadDFloorFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireMixed
)
:
    // The "value" entry is not required
    fvPatchVectorField(p, iF, dict, IOobjectOption::NO_READ),
    refValue_(p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size()),
    source_(p.size(), Zero)
{
    if (!readMixedEntries(dict, requireMixed))
    {
        // Not read (eg, optional and missing): no evaluate possible/need
        return;
    }

    // Could also check/clamp fraction to 0-1 range
    evaluate();
}


// template<class Type>
Foam::breadDFloorFvPatchVectorField::breadDFloorFvPatchVectorField
(
    const breadDFloorFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchVectorField(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper),
    refGrad_(ptf.refGrad_, mapper),
    valueFraction_(ptf.valueFraction_, mapper),
    source_(ptf.source_, mapper)
{
    if (notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


// template<class Type>
Foam::breadDFloorFvPatchVectorField::breadDFloorFvPatchVectorField
(
    const breadDFloorFvPatchVectorField& ptf
)
:
    fvPatchVectorField(ptf),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_),
    source_(ptf.source_)
{}


// template<class Type>
Foam::breadDFloorFvPatchVectorField::breadDFloorFvPatchVectorField
(
    const breadDFloorFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fvPatchVectorField(ptf, iF),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_),
    source_(ptf.source_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class Type>
void Foam::breadDFloorFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchVectorField::autoMap(m);
    refValue_.autoMap(m);
    refGrad_.autoMap(m);
    valueFraction_.autoMap(m);
    source_.autoMap(m);
}


// template<class Type>
void Foam::breadDFloorFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fvPatchVectorField::rmap(ptf, addr);

    const breadDFloorFvPatchVectorField& mptf =
        refCast<const breadDFloorFvPatchVectorField>(ptf);

    refValue_.rmap(mptf.refValue_, addr);
    refGrad_.rmap(mptf.refGrad_, addr);
    valueFraction_.rmap(mptf.valueFraction_, addr);
    source_.rmap(mptf.source_, addr);
}


// template<class Type>
void Foam::breadDFloorFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if(this->db().objectRegistry::foundObject<volScalarField>("impK"))
    {
        const volSymmTensorField& sigma = this->db().objectRegistry::lookupObject<volSymmTensorField>("sigma");
        const volTensorField& gradD = this->db().objectRegistry::lookupObject<volTensorField>("grad(D)");
        const volTensorField& F = this->db().objectRegistry::lookupObject<volTensorField>("F");
        const volScalarField& impK = this->db().objectRegistry::lookupObject<volScalarField>("impK");
        const volScalarField& rImpK = this->db().objectRegistry::lookupObject<volScalarField>("(1|impK)");
        const volVectorField& D = this->db().objectRegistry::lookupObject<volVectorField>("D");
        const surfaceVectorField& Cf = patch().boundaryMesh().mesh().Cf();

        symmTensorField sigmaBound = sigma.boundaryField()[patch().index()];
        tensorField gradDBound = gradD.boundaryField()[patch().index()];
        scalarField impKBound = impK.boundaryField()[patch().index()];
        scalarField rImpKBound = rImpK.boundaryField()[patch().index()];
        tensorField FBound = F.boundaryField()[patch().index()];
        vectorField DBound = D.boundaryField()[patch().index()];
        vectorField CfBound = Cf.boundaryField()[patch().index()];

        tensorField FinvBound = inv(FBound);
        vectorField n = patch().nf();
        vectorField nCurrent = FinvBound.T() & n;
        nCurrent /= mag(nCurrent);

        vectorField gradDForcedBound = (- (nCurrent & sigmaBound) + impKBound * (n & gradDBound)) * rImpKBound;

        forAll(CfBound, faceI)
        {
            if ((CfBound[faceI][0] + DBound[faceI][0]) < floorPos_ )
            {
                this->valueFraction()[faceI] = 1;
                this->refGrad()[faceI] = vector(0,0,0);
                vector oprava = DBound[faceI];
                oprava[0] = floorPos_ - CfBound[faceI][0] * 0.999;
                this->refValue()[faceI] = oprava;
            }
            else
            {
                this->valueFraction()[faceI] = 0;
                this->refGrad()[faceI] = gradDForcedBound[faceI];
            }
        }
    }

    vectorField::operator=
    (
        lerp
        (
            this->patchInternalField() + refGrad_/this->patch().deltaCoeffs(),
            refValue_,
            valueFraction_
        )
    );

    fvPatchVectorField::evaluate();
}


// template<class Type>
Foam::tmp<Foam::vectorField>
Foam::breadDFloorFvPatchVectorField::snGrad() const
{
    return lerp
    (
        refGrad_,
        (refValue_ - this->patchInternalField())*this->patch().deltaCoeffs(),
        valueFraction_
    );
}


// template<class Type>
Foam::tmp<Foam::vectorField>
Foam::breadDFloorFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return vector(pTraits<vector>::one)*(1.0 - valueFraction_);
}


// template<class Type>
Foam::tmp<Foam::vectorField>
Foam::breadDFloorFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return lerp
    (
        refGrad_/this->patch().deltaCoeffs(),
        refValue_,
        valueFraction_
    );
}


// template<class Type>
Foam::tmp<Foam::vectorField>
Foam::breadDFloorFvPatchVectorField::gradientInternalCoeffs() const
{
    return -vector(pTraits<vector>::one)*valueFraction_*this->patch().deltaCoeffs();
}


// template<class Type>
Foam::tmp<Foam::vectorField>
Foam::breadDFloorFvPatchVectorField::gradientBoundaryCoeffs() const
{
    return lerp
    (
        refGrad_,
        this->patch().deltaCoeffs()*refValue_,
        valueFraction_
    );
}

void Foam::breadDFloorFvPatchVectorField::writeScalarEntry(Foam::Ostream& os, Foam::word name, Foam::scalar value) const
{
    os.write("\n\t\t");
    os.write(name);
    os.write("\t");
    os.write(value);
    os.write(";\n");
}


// template<class Type>
void Foam::breadDFloorFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    refValue_.writeEntry("refValue", os);
    refGrad_.writeEntry("refGradient", os);
    valueFraction_.writeEntry("valueFraction", os);
    source_.writeEntry("source", os);
    writeScalarEntry(os, "floorPos", floorPos_);
    fvPatchVectorField::writeValueEntry(os);
}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        breadDFloorFvPatchVectorField
    );
}


// ************************************************************************* //
