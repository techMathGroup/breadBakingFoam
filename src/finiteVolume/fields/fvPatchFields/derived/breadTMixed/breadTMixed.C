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

#include "breadTMixed.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// template<class Type>
bool Foam::breadTMixedFvPatchScalarField::readMixedEntries
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
    dict.readEntry("alpha", alpha_);
    intLamName_ = dict.getOrDefault<word>("intLamName", "lambdaEff");
    TInfDict_ = dict.subDict("TInfTableDict");
    TInfTable_ = interpolationTable<scalar>(TInfDict_);
    this->refValue() = TInfTable_(0);
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// template<class Type>
Foam::breadTMixedFvPatchScalarField::breadTMixedFvPatchScalarField
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
Foam::breadTMixedFvPatchScalarField::breadTMixedFvPatchScalarField
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
Foam::breadTMixedFvPatchScalarField::breadTMixedFvPatchScalarField
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
Foam::breadTMixedFvPatchScalarField::breadTMixedFvPatchScalarField
(
    const breadTMixedFvPatchScalarField& ptf,
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
    alpha_ = ptf.alpha_;
    intLamName_ = ptf.intLamName_;
    TInfTable_ = ptf.TInfTable_;
    TInfDict_ = ptf.TInfDict_;
}


// template<class Type>
Foam::breadTMixedFvPatchScalarField::breadTMixedFvPatchScalarField
(
    const breadTMixedFvPatchScalarField& ptf
)
:
    fvPatchScalarField(ptf),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_),
    source_(ptf.source_)
{
    alpha_ = ptf.alpha_;
    intLamName_ = ptf.intLamName_;
    TInfTable_ = ptf.TInfTable_;
    TInfDict_ = ptf.TInfDict_;
}


// template<class Type>
Foam::breadTMixedFvPatchScalarField::breadTMixedFvPatchScalarField
(
    const breadTMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchScalarField(ptf, iF),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_),
    source_(ptf.source_)
{
    alpha_ = ptf.alpha_;
    intLamName_ = ptf.intLamName_;
    TInfTable_ = ptf.TInfTable_;
    TInfDict_ = ptf.TInfDict_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class Type>
void Foam::breadTMixedFvPatchScalarField::autoMap
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
void Foam::breadTMixedFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fvPatchScalarField::rmap(ptf, addr);

    const breadTMixedFvPatchScalarField& mptf =
        refCast<const breadTMixedFvPatchScalarField>(ptf);

    refValue_.rmap(mptf.refValue_, addr);
    refGrad_.rmap(mptf.refGrad_, addr);
    valueFraction_.rmap(mptf.valueFraction_, addr);
    source_.rmap(mptf.source_, addr);
}


// template<class Type>
void Foam::breadTMixedFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if(this->db().objectRegistry::foundObject<volScalarField>(intLamName_))
    {
        const volScalarField& lambdaEff = this->db().objectRegistry::lookupObject<volScalarField>(intLamName_);
        if (lambdaEff.boundaryField()[this->patch().index()].size() != 0)
        {
            // -- heat transfer to bread computation
            // -- patch deltaCoeffs
            // const volScalarField& lambdaEff = this->db().objectRegistry::lookupObject<volScalarField>(intLamName_);
            // const fvMesh& mesh = patch().boundaryMesh().mesh();
            // const volVectorField& D = this->db().objectRegistry::lookupObject<volVectorField>("D");
            // const surfaceVectorField& Sf = mesh.Sf();

            // vectorField DCells = D.boundaryField()[this->patch().index()].patchInternalField();
            // vectorField DBound = D.boundaryField()[this->patch().index()];
            // vectorField SfBound = Sf.boundaryField()[this->patch().index()];

            // Pout << "DBound size" << DBound.size() <<endl;
            // Pout << "SfBound size" << SfBound.size() <<endl;

            // scalarField Dmag = (DBound - DCells) & SfBound / mag(SfBound);

            // scalarField DCorrect = (DBound - DCells) & mesh
            const scalar t = this->db().time().timeOutputValue();
            scalarField lambdaEffBound = lambdaEff.boundaryField()[this->patch().index()];
            // scalarField f = 1.0 / (1.0 + (lambdaEffBound / (mag(this->patch().delta() + (DBound - DCells)))) / (alpha_));
            scalarField f = 1.0 / (1.0 + (lambdaEffBound * this->patch().deltaCoeffs()) / (alpha_));
            this->valueFraction() = f;
            this->refValue() = TInfTable_(t);
        }
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
Foam::breadTMixedFvPatchScalarField::snGrad() const
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
Foam::breadTMixedFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return scalar(pTraits<scalar>::one)*(1.0 - valueFraction_);
}


// template<class Type>
Foam::tmp<Foam::scalarField>
Foam::breadTMixedFvPatchScalarField::valueBoundaryCoeffs
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
Foam::breadTMixedFvPatchScalarField::gradientInternalCoeffs() const
{
    return -scalar(pTraits<scalar>::one)*valueFraction_*this->patch().deltaCoeffs();
}


// template<class Type>
Foam::tmp<Foam::scalarField>
Foam::breadTMixedFvPatchScalarField::gradientBoundaryCoeffs() const
{
    return lerp
    (
        refGrad_,
        this->patch().deltaCoeffs()*refValue_,
        valueFraction_
    );
}

void Foam::breadTMixedFvPatchScalarField::writeScalarEntry(Foam::Ostream& os, Foam::word name, Foam::scalar value) const
{
    os.write("\n\t\t");
    os.write(name);
    os.write("\t");
    os.write(value);
    os.write(";\n");
}

// template<class Type>
void Foam::breadTMixedFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    refValue_.writeEntry("refValue", os);
    refGrad_.writeEntry("refGradient", os);
    valueFraction_.writeEntry("valueFraction", os);
    breadTMixedFvPatchScalarField::writeScalarEntry(os, "alpha", alpha_);
    source_.writeEntry("source", os);
    os.writeEntry("TInfTableDict", TInfDict_);
    fvPatchScalarField::writeValueEntry(os);
}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        breadTMixedFvPatchScalarField
    );
}


// ************************************************************************* //
