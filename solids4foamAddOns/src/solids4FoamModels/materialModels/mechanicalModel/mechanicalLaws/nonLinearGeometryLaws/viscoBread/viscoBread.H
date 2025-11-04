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
    viscoBread

Description
    Hyper-elasto-plastic mechanical law as described by Simo & Hughes (1998)
    in Box 9.1.

    Neo-Hookean elasticity with Mises/J2 plasticity.

    For elastic parameters, the user can specify:
        - Young's modulus (E) and Poisson's ratio (nu)
    or
        - Shear modulus (mu) and bulk modulus (K)

    More details found in:

    Simo & Hughes, Computational Inelasticity, 1998, Springer.

    P. Cardiff, Z. Tukovic, P. De Jaeger, M. Clancy and A. Ivankovic. A
    Lagrangian cell-centred finite volume method for metal forming simulation,
    doi=10.1002/nme.5345.

SourceFiles
    viscoBread.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef viscoBread_H
#define viscoBread_H

#include "mechanicalLaw.H"
#include "interpolationTable.H"
#ifdef OPENFOAM_NOT_EXTEND
    #include "surfaceFields.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class viscoBread
:
    public mechanicalLaw
{
    // Private data

        //- Shear modulus
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar K_;

        //- Jacobian of the toal deformation gradient volField
        autoPtr<volScalarField> JPtr_;

        //- Jacobian of the toal deformation gradient surfaceField
        autoPtr<surfaceScalarField> JfPtr_;

        //- Total plastic strain
        volSymmTensorField epsilonP_;

        //- Total plastic strain surface field
        surfaceSymmTensorField epsilonPf_;

        //- Incremental change of plastic strain
        volSymmTensorField DEpsilonP_;

        //- Incremental change of plastic strain surface field
        surfaceSymmTensorField DEpsilonPf_;

        //- Elastic left Cauchy-Green trial strain tensor with volumetric term
        //  removed. Trial means that plasticity is neglected.
        volSymmTensorField bEbarTrial_;

        //- Elastic left Cauchy-Green trial strain tensor with volumetric term
        //  removed surface field
        surfaceSymmTensorField bEbarTrialf_;

        //- Elastic left Cauchy-Green strain tensor with volumetric term removed
        volSymmTensorField bEbar_;

        //- Elastic left Cauchy-Green strain tensor with volumetric term removed
        //  surface field
        surfaceSymmTensorField bEbarf_;

        //- Plastic multiplier increment - plastric strain scaled by sqrt(2/3)
        volScalarField DLambda_;

        //- Plastic multiplier increment - surface field
        surfaceScalarField DLambdaf_;

        //- plasticN is the return direction to the yield surface
        volSymmTensorField plasticN_;

        //- plasticN surface field
        surfaceSymmTensorField plasticNf_;

        //- Maximum allowed error in the plastic strain integration
        const scalar maxDeltaErr_;

        //- Tolerance for Newton loop
        static scalar LoopTol_;

        //- Maximum number of iterations for Newton loop
        static label MaxNewtonIter_;

        //- finiteDiff is the delta for finite difference differentiation
        static scalar finiteDiff_;

        //- Store sqrt(2/3) as it is used often
        static scalar sqrtTwoOverThree_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        viscoBread(const viscoBread&);

        //- Disallow default bitwise assignment
        void operator=(const viscoBread&);

        //- Make the J field
        void makeJ();

        //- Return a reference to the J field
        volScalarField& J();

        //- Make the Jf field
        void makeJf();

        //- Return a reference to the Jf field
        surfaceScalarField& Jf();

public:

    //- Runtime type information
    TypeName("viscoBread");

    // Static data members


    // Constructors

        //- Construct from dictionary
        viscoBread
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~viscoBread();


    // Member Functions

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Update the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Update the stress surface field
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update the yield stress: called at end of time-step
        virtual void updateTotalFields();

        //- Return the desired new time-step
        virtual scalar newDeltaT();

        //- Set restart: some mechanical laws may change the writeOpt of fields
        //  accordingly
        virtual void setRestart();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
