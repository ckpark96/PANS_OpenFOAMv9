/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "PANSkEpsilon.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void PANSkEpsilon<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(kU_)/epsilonU_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
    // BasicMomentumTransportModel::correctNut();
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> PANSkEpsilon<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            kU_,
            dimVolume*this->rho_.dimensions()*kU_.dimensions()
            /dimTime
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> PANSkEpsilon<BasicMomentumTransportModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilonU_,
            dimVolume*this->rho_.dimensions()*epsilonU_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
PANSkEpsilon<BasicMomentumTransportModel>::PANSkEpsilon
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),
    // start of new stuff
    fEpsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fEpsilon",
            this->coeffDict_,
            1.0
        )
    ),
    // Fixed fK
    fK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fK",
            this->coeffDict_,
            1.0 // to be varied based on preliminary RANS calculation
        )
    ),
    // // Variable fK
    // uLim_
    // (
    //     dimensioned<scalar>::lookupOrAddToDict
    //     (
    //         "fKupperLimit",
    //         this->coeffDict_,
    //         1.0
    //     )
    // ),
    // loLim_
    // (
    //     dimensioned<scalar>::lookupOrAddToDict
    //     (
    //         "fKlowerLimit",
    //         this->coeffDict_,
    //         0.1
    //     )
    // ),
    // fK_
    // (
    //     IOobject
    //     (
    //         IOobject::groupName("fK", alphaRhoPhi.group()),
    //         this->runTime_.timeName(),
    //         this->mesh_,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     this->mesh_,
    //     dimensionedScalar("zero",loLim_)
    // ),

    C2U_
    (C1_ + (fK_/fEpsilon_)*(C2_ - C1_)),
        
        // IOobject
        // (
        //     "C2U",
        //     this->runTime_.timeName(),
        //     this->mesh_,
        //     IOobject::NO_READ,
        //     IOobject::NO_WRITE
        // ),
        
    // delta_
    // (
    //     LESdelta::New
    //     (
    //         IOobject::groupName("delta", U.group()),
    //         *this,
    //         this->coeffDict_
    //     )
    // ),
// end of new stuff
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    kU_
    (
        IOobject
        (
            IOobject::groupName("kU", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        k_*fK_,
        k_.boundaryField().types()
    ),
    epsilonU_
    (
        IOobject
        (
            IOobject::groupName("epsilonU", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
    	),
    	epsilon_*fEpsilon_,
        epsilon_.boundaryField().types()
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);
    bound(kU_, fK_*this->kMin_);
    bound(epsilonU_, fEpsilon_*this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
        correctNut();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool PANSkEpsilon<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        fEpsilon_.readIfPresent(this->coeffDict());
        // uLim_.readIfPresent(this->coeffDict());
        // loLim_.readIfPresent(this->coeffDict());
        fK_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void PANSkEpsilon<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField::Internal G
    (
        this->GName(),
        nut()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );
    tgradU.clear();

    // Update epsilon and G at the wall
    epsilonU_.boundaryFieldRef().updateCoeffs();


    // Dissipation equation
    tmp<fvScalarMatrix> epsUEqn
    (
        fvm::ddt(alpha, rho, epsilonU_)
      + fvm::div(alphaRhoPhi, epsilonU_)
      - fvm::laplacian(alpha*rho*DepsilonUEff(), epsilonU_)
     ==
        C1_*alpha()*rho()*G*epsilonU_()/kU_()
      - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha()*rho()*divU, epsilonU_)
      - fvm::Sp(C2U_*alpha()*rho()*epsilonU_()/kU_(), epsilonU_)
      + epsilonSource()
      + fvModels.source(alpha, rho, epsilonU_)
    );

    epsUEqn.ref().relax();
    fvConstraints.constrain(epsUEqn.ref());
    epsUEqn.ref().boundaryManipulate(epsilonU_.boundaryFieldRef());
    solve(epsUEqn);
    fvConstraints.constrain(epsilonU_);
    bound(epsilonU_, fEpsilon_*this->epsilonMin_);

    // Unresolved turbulent kinetic energy equation
    tmp<fvScalarMatrix> kUEqn
    (
        fvm::ddt(alpha, rho, kU_)
      + fvm::div(alphaRhoPhi, kU_)
      - fvm::laplacian(alpha*rho*DkUEff(), kU_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, kU_)
      - fvm::Sp(alpha()*rho()*epsilonU_()/kU_(), kU_)
      + kSource()
      + fvModels.source(alpha, rho, kU_)
    );

    kUEqn.ref().relax();
    fvConstraints.constrain(kUEqn.ref());
    solve(kUEqn);
    fvConstraints.constrain(kU_);
    bound(kU_, fK_*this->kMin_);

    // Calculation of Turbulent kinetic energy and Dissipation rate
    k_ = kU_/fK_;
    k_.correctBoundaryConditions();
    
    epsilon_ = epsilonU_/fEpsilon_;
    epsilon_.correctBoundaryConditions();

    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    correctNut();


// For Variable fK
    // Recalculate fK and C2U with new kU and epsilonU

    // // Calculate the turbulence integral length scale
    // volScalarField::Internal Lambda
    // (
    // 	pow(k_,1.5)/epsilon_
    // );

    // // update fK
    // fK_.primitiveFieldRef() = min(max(
    // 	sqrt(Cmu_.value())*pow(delta_()/Lambda,2.0/3.0), loLim_.value()), uLim_.value());
    
    // // update C2U
    // C2U = C1_+ (fK_/fEpsilon_)*(C2_-C1_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
