/*---------------------------------------------------------------------------*\
SourceFiles
    PANSkEpsilonPANS.C

\*---------------------------------------------------------------------------*/

#ifndef PANSkEpsilon_H
#define PANSkEpsilon_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class PANSkEpsilon Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class PANSkEpsilon
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel>>
{
protected:

    // Protected data

        // Model coefficients
            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar fEpsilon_;
            //dimensionedScalar uLim_;
            //dimensionedScalar loLim_;
            dimensionedScalar fK_;
            dimensionedScalar C2U_;

        // Fields
            // volScalarField fK_;
            // volScalarField C2U;
            //- Run-time selectable delta model
            // autoPtr<Foam::LESdelta> delta_;
            volScalarField k_;
            volScalarField epsilon_;
            volScalarField kU_;
            volScalarField epsilonU_;
            


    // Protected Member Functions

        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("PANSkEpsilon");


    // Constructors

        //- Construct from components
        PANSkEpsilon
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            //const word& propertiesName = momentumTransportModel::propertiesName,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        PANSkEpsilon(const PANSkEpsilon&) = delete;


    //- Destructor
    virtual ~PANSkEpsilon()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                (this->nut_/sigmak_ + this->nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return volScalarField::New
            (
                "DepsilonEff",
                (this->nut_/sigmaEps_ + this->nu())
            );
        }

        //- Return the effective diffusivity for unresolved k
        tmp<volScalarField> DkUEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                ("DkUEff", this->nut_/(fK_*fK_*sigmak_/fEpsilon_)
                 + this->nu())
            );
        }
        //- Return the effective diffusivity for unresolved epsilon
        tmp<volScalarField> DepsilonUEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                ("DepsilonUEff", this->nut_/(fK_*fK_*sigmaEps_/fEpsilon_)
                 + this->nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();



    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const PANSkEpsilon&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PANSkEpsilon.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
