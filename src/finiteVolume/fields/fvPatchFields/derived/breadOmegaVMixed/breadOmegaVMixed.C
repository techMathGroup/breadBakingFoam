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

#include "breadOmegaVMixed.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// template<class Type>
bool Foam::breadOmegaVMixedFvPatchScalarField::readMixedEntries
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
    omegaVInfDict_ = dict.subDict("omegaVInfTableDict");
    omegaVInfTable_ = interpolationTable<scalar>(omegaVInfDict_);
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// template<class Type>
Foam::breadOmegaVMixedFvPatchScalarField::breadOmegaVMixedFvPatchScalarField
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
Foam::breadOmegaVMixedFvPatchScalarField::breadOmegaVMixedFvPatchScalarField
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
Foam::breadOmegaVMixedFvPatchScalarField::breadOmegaVMixedFvPatchScalarField
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
Foam::breadOmegaVMixedFvPatchScalarField::breadOmegaVMixedFvPatchScalarField
(
    const breadOmegaVMixedFvPatchScalarField& ptf,
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
    kM_ = ptf.kM_;
    omegaVInfTable_ = ptf.omegaVInfTable_;
    omegaVInfDict_ = ptf.omegaVInfDict_;
}


// template<class Type>
Foam::breadOmegaVMixedFvPatchScalarField::breadOmegaVMixedFvPatchScalarField
(
    const breadOmegaVMixedFvPatchScalarField& ptf
)
:
    fvPatchScalarField(ptf),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_),
    source_(ptf.source_)
{
    kM_ = ptf.kM_;
    omegaVInfTable_ = ptf.omegaVInfTable_;
    omegaVInfDict_ = ptf.omegaVInfDict_;
}


// template<class Type>
Foam::breadOmegaVMixedFvPatchScalarField::breadOmegaVMixedFvPatchScalarField
(
    const breadOmegaVMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchScalarField(ptf, iF),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_),
    source_(ptf.source_)
{
    kM_ = ptf.kM_;
    omegaVInfTable_ = ptf.omegaVInfTable_;
    omegaVInfDict_ = ptf.omegaVInfDict_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class Type>
void Foam::breadOmegaVMixedFvPatchScalarField::autoMap
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
void Foam::breadOmegaVMixedFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fvPatchScalarField::rmap(ptf, addr);

    const breadOmegaVMixedFvPatchScalarField& mptf =
        refCast<const breadOmegaVMixedFvPatchScalarField>(ptf);

    refValue_.rmap(mptf.refValue_, addr);
    refGrad_.rmap(mptf.refGrad_, addr);
    valueFraction_.rmap(mptf.valueFraction_, addr);
    source_.rmap(mptf.source_, addr);
}


// template<class Type>
void Foam::breadOmegaVMixedFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    if(this->db().objectRegistry::foundObject<volScalarField>("DEff"))
    {
        const volScalarField& DEff = this->db().objectRegistry::lookupObject<volScalarField>("DEff");
        if (DEff.boundaryField()[this->patch().index()].size() != 0)
        {
            const scalar t = this->db().time().timeOutputValue();
            // const fvMesh& mesh = patch().boundaryMesh().mesh();
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
            const volScalarField& Mg = this->db().objectRegistry::lookupObject<volScalarField>("Mg");
            const volScalarField& pG = this->db().objectRegistry::lookupObject<volScalarField>("pG");
            const volScalarField& alphaS = this->db().objectRegistry::lookupObject<volScalarField>("alphaS");
            const volScalarField& alphaL = this->db().objectRegistry::lookupObject<volScalarField>("alphaL");
            const volScalarField& rhoG = this->db().objectRegistry::lookupObject<volScalarField>("rhoG");
            // const volVectorField& D = this->db().objectRegistry::lookupObject<volVectorField>("D");
            // const surfaceVectorField& Sf = mesh.Sf();

            // vectorField DCells = D.boundaryField()[this->patch().index()].patchInternalField();
            // vectorField DBound = D.boundaryField()[this->patch().index()];
            // vectorField SfBound = Sf.boundaryField()[this->patch().index()];
            // scalarField Dmag =  (DBound - DCells) & SfBound / mag(SfBound);

            // Pout << "min Dmag" <<min(Dmag) << "max(Dmag)" << max(Dmag) <<endl;

            scalarField rhoGBound = rhoG.boundaryField()[this->patch().index()];
            scalarField permBound = perm.boundaryField()[this->patch().index()];
            scalarField DEffmBound = DEff.boundaryField()[this->patch().index()];
            scalarField MgBound = Mg.boundaryField()[this->patch().index()];
            scalarField MgCells = Mg.boundaryField()[this->patch().index()].patchInternalField();
            scalarField TBound = T.boundaryField()[this->patch().index()];
            scalarField pGBound = pG.boundaryField()[this->patch().index()];
            scalarField pGCells = pG.boundaryField()[this->patch().index()].patchInternalField();
                    
            scalarField alphaGBound = 1 - alphaS.boundaryField()[this->patch().index()] - alphaL.boundaryField()[this->patch().index()];
            scalarField K1Bound = rhoGBound * permBound * this->patch().deltaCoeffs();
            scalarField K2Bound = rhoGBound * DEff * this->patch().deltaCoeffs();
            // scalarField K1Bound = rhoGBound * permBound / (mag(this->patch().delta() + (DBound - DCells)));
            // scalarField K2Bound = rhoGBound * DEff / (mag(this->patch().delta() + (DBound - DCells)));

            // scalarField denominator = K1Bound * (pGBound - pGBound) + K2Bound + K2Bound / MgBound * (MgBound - MgCells) + kM_ * rhoGBound * alphaGBound;
            // scalarField denominator = K1Bound * (pGBound - pGBound) + K2Bound + K2Bound / MgBound * (MgBound - MgCells) + kM_ * rhoGBound;
            scalarField denominator = K1Bound * (pGBound - pGCells) + K2Bound + kM_ * rhoGBound;
            scalarField f = 1.0 - K2Bound / denominator;
            scalarField a = (kM_ * rhoGBound * omegaVInfTable_(t)) / denominator;

            this->valueFraction() = f;
            this->refValue() = a / f;
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
Foam::breadOmegaVMixedFvPatchScalarField::snGrad() const
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
Foam::breadOmegaVMixedFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return scalar(pTraits<scalar>::one)*(1.0 - valueFraction_);
}


// template<class Type>
Foam::tmp<Foam::scalarField>
Foam::breadOmegaVMixedFvPatchScalarField::valueBoundaryCoeffs
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
Foam::breadOmegaVMixedFvPatchScalarField::gradientInternalCoeffs() const
{
    return -scalar(pTraits<scalar>::one)*valueFraction_*this->patch().deltaCoeffs();
}


// template<class Type>
Foam::tmp<Foam::scalarField>
Foam::breadOmegaVMixedFvPatchScalarField::gradientBoundaryCoeffs() const
{
    return lerp
    (
        refGrad_,
        this->patch().deltaCoeffs()*refValue_,
        valueFraction_
    );
}

void Foam::breadOmegaVMixedFvPatchScalarField::writeScalarEntry(Foam::Ostream& os, Foam::word name, Foam::scalar value) const
{
    os.write("\n\t\t");
    os.write(name);
    os.write("\t");
    os.write(value);
    os.write(";\n");
}


// template<class Type>
void Foam::breadOmegaVMixedFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    refValue_.writeEntry("refValue", os);
    refGrad_.writeEntry("refGradient", os);
    valueFraction_.writeEntry("valueFraction", os);
    source_.writeEntry("source", os);
    os.writeEntry("kM", this->kM_);
    os.writeEntry("omegaVInfTableDict", omegaVInfDict_);
    // writeScalarEntry(os, "kM", kM_);
    // writeScalarEntry(os, "omegaVInf", omegaVInf_);
    fvPatchScalarField::writeValueEntry(os);
}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        breadOmegaVMixedFvPatchScalarField
    );
}


// ************************************************************************* //
