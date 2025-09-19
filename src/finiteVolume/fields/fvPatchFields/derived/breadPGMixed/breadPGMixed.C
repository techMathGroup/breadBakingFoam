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

#include "breadPGMixed.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// template<class Type>
bool Foam::breadPGMixedFvPatchScalarField::readMixedEntries
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
    dict.readEntry("kM", kM_);
    dict.readEntry("pGInf", pGInf_);
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// template<class Type>
Foam::breadPGMixedFvPatchScalarField::breadPGMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchScalarField(p, iF),
    refValue_(p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size()),
    source_(p.size(), Zero)
{}


// template<class Type>
Foam::breadPGMixedFvPatchScalarField::breadPGMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const Foam::zero
)
:
    fvPatchScalarField(p, iF),
    refValue_(p.size(), Zero),
    refGrad_(p.size(), Zero),
    valueFraction_(p.size(), Zero),
    source_(p.size(), Zero)
{}


// template<class Type>
Foam::breadPGMixedFvPatchScalarField::breadPGMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireMixed
)
:
    // The "value" entry is not required
    fvPatchScalarField(p, iF, dict, IOobjectOption::NO_READ),
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
Foam::breadPGMixedFvPatchScalarField::breadPGMixedFvPatchScalarField
(
    const breadPGMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchScalarField(ptf, p, iF, mapper),
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
Foam::breadPGMixedFvPatchScalarField::breadPGMixedFvPatchScalarField
(
    const breadPGMixedFvPatchScalarField& ptf
)
:
    fvPatchScalarField(ptf),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_),
    source_(ptf.source_)
{}


// template<class Type>
Foam::breadPGMixedFvPatchScalarField::breadPGMixedFvPatchScalarField
(
    const breadPGMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchScalarField(ptf, iF),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_),
    source_(ptf.source_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class Type>
void Foam::breadPGMixedFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchScalarField::autoMap(m);
    refValue_.autoMap(m);
    refGrad_.autoMap(m);
    valueFraction_.autoMap(m);
    source_.autoMap(m);
}


// template<class Type>
void Foam::breadPGMixedFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fvPatchScalarField::rmap(ptf, addr);

    const breadPGMixedFvPatchScalarField& mptf =
        refCast<const breadPGMixedFvPatchScalarField>(ptf);

    refValue_.rmap(mptf.refValue_, addr);
    refGrad_.rmap(mptf.refGrad_, addr);
    valueFraction_.rmap(mptf.valueFraction_, addr);
    source_.rmap(mptf.source_, addr);
}


// template<class Type>
void Foam::breadPGMixedFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if(this->db().objectRegistry::foundObject<volScalarField>("permGLViscG"))
    {
        IOdictionary transportProperties = this->db().objectRegistry::lookupObject<IOdictionary>("transportProperties");
        IOdictionary thermophysicalProperties = this->db().objectRegistry::lookupObject<IOdictionary>("thermophysicalProperties");
        // -- heat transfer to bread computation
        // -- patch deltaCoeffs
        scalar molMRef;
        // dimensionedScalar permGLViscG, univR;
        dimensionedScalar univR;
        // transportProperties.readEntry("kGOver",kG);
        transportProperties.subDict("genProps").readEntry("univR",univR);

        // thermophysicalProperties.subDict("solid").readEntry("permGLViscG",permGLViscG);
        thermophysicalProperties.subDict("mixture").subDict("specie").readEntry("molWeight",molMRef);
        molMRef = molMRef * 1e-3;

        const volScalarField& perm = this->db().objectRegistry::lookupObject<volScalarField>("permGLViscG");
        const volScalarField& T = this->db().objectRegistry::lookupObject<volScalarField>("T");
        const volScalarField& rhoG = this->db().objectRegistry::lookupObject<volScalarField>("rhoG");
        const volScalarField& Mg = this->db().objectRegistry::lookupObject<volScalarField>("Mg");
        const volScalarField& alphaS = this->db().objectRegistry::lookupObject<volScalarField>("alphaS");
        const volScalarField& alphaL = this->db().objectRegistry::lookupObject<volScalarField>("alphaL");

        scalarField rhoGBound = rhoG.boundaryField()[patch().index()];
        scalarField permBound = perm.boundaryField()[patch().index()];
        scalarField MgBound = Mg.boundaryField()[patch().index()];
        scalarField TBound = T.boundaryField()[patch().index()];
        
        scalarField alphaGBound = 1 - alphaS.boundaryField()[patch().index()] - alphaL.boundaryField()[patch().index()];
        scalarField K1Bound = rhoGBound * permBound * patch().deltaCoeffs();
        scalarField TInf = TBound;

        // scalarField denominator = K1Bound + kM_ * MgBound / univR.value() / TBound * alphaGBound;
        scalarField denominator = K1Bound + kM_ * MgBound / univR.value() / TBound;
        scalarField f = 1 - K1Bound / denominator;
        scalarField a = (kM_ * molMRef * pGInf_ / (univR.value() * TInf)) / denominator;

        this->valueFraction() = f;
        this->refValue() = a / f;
    }

    scalarField::operator=
    (
        lerp
        (
            this->patchInternalField() + refGrad_/this->patch().deltaCoeffs(),
            refValue_,
            valueFraction_
        )
    );

    fvPatchScalarField::evaluate();
}


// template<class Type>
Foam::tmp<Foam::scalarField>
Foam::breadPGMixedFvPatchScalarField::snGrad() const
{
    return lerp
    (
        refGrad_,
        (refValue_ - this->patchInternalField())*this->patch().deltaCoeffs(),
        valueFraction_
    );
}


// template<class Type>
Foam::tmp<Foam::scalarField>
Foam::breadPGMixedFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return scalar(pTraits<scalar>::one)*(1.0 - valueFraction_);
}


// template<class Type>
Foam::tmp<Foam::scalarField>
Foam::breadPGMixedFvPatchScalarField::valueBoundaryCoeffs
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
Foam::tmp<Foam::scalarField>
Foam::breadPGMixedFvPatchScalarField::gradientInternalCoeffs() const
{
    return -scalar(pTraits<scalar>::one)*valueFraction_*this->patch().deltaCoeffs();
}


// template<class Type>
Foam::tmp<Foam::scalarField>
Foam::breadPGMixedFvPatchScalarField::gradientBoundaryCoeffs() const
{
    return lerp
    (
        refGrad_,
        this->patch().deltaCoeffs()*refValue_,
        valueFraction_
    );
}

void Foam::breadPGMixedFvPatchScalarField::writeScalarEntry(Foam::Ostream& os, Foam::word name, Foam::scalar value) const
{
    os.write("\n\t\t");
    os.write(name);
    os.write("\t");
    os.write(value);
    os.write(";\n");
}

// template<class Type>
void Foam::breadPGMixedFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    refValue_.writeEntry("refValue", os);
    refGrad_.writeEntry("refGradient", os);
    valueFraction_.writeEntry("valueFraction", os);
    source_.writeEntry("source", os);
    writeScalarEntry(os, "kM", kM_);
    writeScalarEntry(os, "pGInf", pGInf_);
    fvPatchScalarField::writeValueEntry(os);
}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        breadPGMixedFvPatchScalarField
    );
}


// ************************************************************************* //
