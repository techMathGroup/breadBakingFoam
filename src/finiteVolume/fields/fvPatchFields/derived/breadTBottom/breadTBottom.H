/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::breadTBottomFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a base class for 'breadTBottom' type boundary
    conditions, i.e. conditions that mix fixed value and patch-normal gradient
    conditions.

    Usage:

    type breadTBottom;
    alpha           10; // -- heat transfer coefficient
    TInf            463;    // -- oven temperature

    // -- mixed BC mandatory entires
    refValue        uniform 303;
    refGradient     uniform 0;
    valueFraction   uniform 0;
    value           uniform 303;

Note
    Boundary condition for temperature field for bread in oven.

SourceFiles
    breadTBottomFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef breadTBottomFvPatchScalarField_H
#define breadTBottomFvPatchScalarField_H

#include "fvPatchField.H"
#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "updateableSnGrad.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class breadTBottomFvPatchField Declaration
\*---------------------------------------------------------------------------*/

// template<class Type>
class breadTBottomFvPatchScalarField
:
    public fvPatchScalarField
{
    // Private Data

        //- Value field
        scalarField refValue_;

        //- Normal gradient field
        scalarField refGrad_;

        //- Fraction (0-1) of value used for boundary condition
        scalarField valueFraction_;

        //- Source field
        scalarField source_;

        //- heat transfer coefficient
        scalar alpha_;

        //- internal heat transfer coefficient field name 
        word intLamName_;

        //- dictionary with T in oven
        dictionary TInfDict_;

        //- Interpolation table for temperature vs time
        interpolationTable<scalar> qInfTable_;

protected:

        //- Read the 'refValue', 'refGradient' and 'valueFraction' entries
        //- into their respective places.
        //  The reading can be optional (default), mandatory etc.
        //  If 'refValue' is to be read, refGradient and valueFraction must
        //  also exist.
        //  \returns True on success
        bool readMixedEntries
        (
            const dictionary& dict,
            IOobjectOption::readOption readOpt = IOobjectOption::LAZY_READ
        );


public:

    //- Runtime type information
    TypeName("breadTBottom");


    // Constructors

        //- Construct from patch and internal field
        breadTBottomFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch and internal field,
        //- initialise as zero-gradient
        breadTBottomFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const Foam::zero
        );

        //- Construct from patch, internal field and dictionary
        breadTBottomFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            //! The "refValue", "refGradient", "valueFraction" entries
            //! (default: mandatory)
            IOobjectOption::readOption requireMixed = IOobjectOption::MUST_READ
        );

        //- Construct by mapping the given breadTBottomFvPatchField onto a new patch
        breadTBottomFvPatchScalarField
        (
            const breadTBottomFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        breadTBottomFvPatchScalarField
        (
            const breadTBottomFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new breadTBottomFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        breadTBottomFvPatchScalarField
        (
            const breadTBottomFvPatchScalarField&,
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
                new breadTBottomFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- True: this patch field fixes a value.
        virtual bool fixesValue() const { return true; }

        //- False: this patch field is not altered by assignment.
        virtual bool assignable() const { return false; }


        // Return defining fields

            virtual scalarField& refValue()
            {
                return refValue_;
            }

            virtual const scalarField& refValue() const
            {
                return refValue_;
            }

            virtual scalarField& refGrad()
            {
                return refGrad_;
            }

            virtual const scalarField& refGrad() const
            {
                return refGrad_;
            }

            virtual scalarField& valueFraction()
            {
                return valueFraction_;
            }

            virtual const scalarField& valueFraction() const
            {
                return valueFraction_;
            }

            virtual scalarField& source()
            {
                return source_;
            }

            virtual const scalarField& source() const
            {
                return source_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<scalarField> snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<scalarField> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<scalarField> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<scalarField> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<scalarField> gradientBoundaryCoeffs() const;


        //- Write
        virtual void write(Ostream&) const;

        void writeScalarEntry(Foam::Ostream&, Foam::word, Foam::scalar) const;


    // Member operators

        virtual void operator=(const UList<scalar>&) {}

        virtual void operator=(const fvPatchScalarField&) {}
        virtual void operator+=(const fvPatchScalarField&) {}
        virtual void operator-=(const fvPatchScalarField&) {}
        virtual void operator*=(const fvPatchScalarField&) {}
        virtual void operator/=(const fvPatchScalarField&) {}

        virtual void operator+=(const scalarField&) {}
        virtual void operator-=(const scalarField&) {}

        virtual void operator*=(const scalarField&) {}
        virtual void operator/=(const scalarField&) {}

        virtual void operator=(const scalar&) {}
        virtual void operator+=(const scalar&) {}
        virtual void operator-=(const scalar&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "breadTBottomFvPatchField.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
