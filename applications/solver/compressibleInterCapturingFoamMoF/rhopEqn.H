{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", rhof*fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        rhof*fvc::flux(HbyA)
      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, rhoPhi))
    );
    MRF.makeRelative(rhof, phiHbyA);

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rhorAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

    tmp<fvScalarMatrix> p_rghEqnComp1;
    tmp<fvScalarMatrix> p_rghEqnComp2;

    // if (pimple.transonic())
    // {
    //     #include "rhofs.H"

    //     surfaceScalarField phid1("phid1", fvc::interpolate(psi1)*phi);
    //     surfaceScalarField phid2("phid2", fvc::interpolate(psi2)*phi);

    //     p_rghEqnComp1 =
    //         pos(alpha1)
    //        *(
    //             (
    //                 fvc::ddt(alpha1, rho1) + fvc::div(alphaPhi1*rho1f)
    //               - (fvOptions(alpha1, mixture.thermo1().rho())&rho1)
    //             )/rho1
    //           - fvc::ddt(alpha1) - fvc::div(alphaPhi1)
    //           + (alpha1/rho1)
    //            *correction
    //             (
    //                 psi1*fvm::ddt(p_rgh)
    //               + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
    //             )
    //         );
    //     p_rghEqnComp1.ref().relax();

    //     p_rghEqnComp2 =
    //         pos(alpha2)
    //        *(
    //            (
    //                fvc::ddt(alpha2, rho2) + fvc::div(alphaPhi2*rho2f)
    //              - (fvOptions(alpha2, mixture.thermo2().rho())&rho2)
    //            )/rho2
    //          - fvc::ddt(alpha2) - fvc::div(alphaPhi2)
    //          + (alpha2/rho2)
    //           *correction
    //            (
    //                psi2*fvm::ddt(p_rgh)
    //              + fvm::div(phid2, p_rgh) - fvm::Sp(fvc::div(phid2), p_rgh)
    //            )
    //        );
    //     p_rghEqnComp2.ref().relax();
    // }
    // else
    {
        #include "rhofs.H"

        p_rghEqnComp1 =
            alpha1*(fvc::ddt(rho1) + psi1*correction(fvm::ddt(p_rgh)));

        p_rghEqnComp2 =
            alpha2*(fvc::ddt(rho2) + psi2*correction(fvm::ddt(p_rgh)));
    }

    // Cache p_rgh prior to solve for density update
    volScalarField p_rgh_0(p_rgh);

    // Info << "p_rghEqnComp2.dimensions " << p_rghEqnComp2.dimensions() << endl;

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            p_rghEqnComp1 + p_rghEqnComp2
          + fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, p_rgh)
        );

        solve
        (
            pEqn,
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            p = max(p_rgh + (alpha1*rho1 + alpha2*rho2)*gh, pMin);
            p_rgh = p - (alpha1*rho1 + alpha2*rho2)*gh;


            rhoPhi = phiHbyA + pEqn.flux();
            phi = rhoPhi/rhof;

            U = HbyA+ rAU*fvc::reconstruct((phig + pEqn.flux())/rhorAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    // Update densities from change in p_rgh
    mixture.thermo1().correctRho(psi1*(p_rgh - p_rgh_0));
    mixture.thermo2().correctRho(psi2*(p_rgh - p_rgh_0));

    volVectorField rhoPhi2("rhoPhi2",fvc::reconstruct(rhoPhi));
    volVectorField phi2("phi2",fvc::reconstruct(phi));

    rhoPhi2.write();
    phi2.write();

    // rho = alpha1*rho1 + alpha2*rho2;
    {
       #include "compressibleAlphaRhoEqnSubCycle.H"
    }

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();

    K = 0.5*magSqr(U);
}
