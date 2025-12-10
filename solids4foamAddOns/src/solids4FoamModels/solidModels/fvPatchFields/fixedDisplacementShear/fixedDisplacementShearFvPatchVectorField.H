/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    fixedDisplacementShearFvPatchVectorField

Description
    The component of the displacement normal to the patch is applied,
    and the shear traction is set to zero.

    Patch contains non-orthogonal correction.

    Optionally time-varying by specifying a time series.

SourceFiles
    fixedDisplacementShearFvPatchVectorField.C

Author
    Philip Cardiff, UCD/UT. All rights reserved.
    Based on velocityDisplacement by Aleksandar Karac

\*---------------------------------------------------------------------------*/

#ifndef fixedDisplacementShearFvPatchVectorField_H
#define fixedDisplacementShearFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidDirectionMixedFvPatchVectorField.H"
#include "Switch.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedDisplacementShearFvPatchVectorField
:
    public solidDirectionMixedFvPatchVectorField
{
    // Private data

        //- Constant total displacement field
        vectorField totalDisp_;

        //- Displacement time series
        interpolationTable<vector> dispSeries_;

        //- Switch to force the shear gradient to be zero
        const Switch forceZeroShearGrad_;

public:

    //- Runtime type information
    TypeName("fixedDisplacementShear");


    // Constructors

        //- Construct from patch and internal field
        fixedDisplacementShearFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedDisplacementShearFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // fixedDisplacementShearFvPatchVectorField
        //  onto a new patch
        fixedDisplacementShearFvPatchVectorField
        (
            const fixedDisplacementShearFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

#ifdef OPENFOAM_ORG
        //- Construct as copy
        fixedDisplacementShearFvPatchVectorField
        (
            const fixedDisplacementShearFvPatchVectorField&
        ) = delete;
#else
        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new fixedDisplacementShearFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        fixedDisplacementShearFvPatchVectorField
        (
            const fixedDisplacementShearFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new fixedDisplacementShearFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~fixedDisplacementShearFvPatchVectorField()
        {}


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
