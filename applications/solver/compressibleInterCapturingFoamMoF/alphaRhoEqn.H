word alphaRhoScheme("div(phi,alphaRho)");
word alpharRhoScheme("div(phirb,alphaRho)");

// Standard face-flux compression coefficient
surfaceScalarField rhoPhi1(phi*rho1f);
surfaceScalarField rhoPhic(mixture.cAlpha()*mag(rhoPhi1/mesh.magSf()));


for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
{

    surfaceScalarField rhoPhir(rhoPhic*mixture.nHatf());

    alphaRhoPhi1 =
        fvc::flux
        (
            rhoPhi1,
            alpha1,
            alphaRhoScheme
        )
        + fvc::flux
        (
            -fvc::flux(-rhoPhir, alpha2, alpharRhoScheme),
            alpha1,
            alpharRhoScheme
        );

    // if (MULESCorr)
    // {
    //     tmp<surfaceScalarField> talphaPhi1Corr(alphaPhiUn - alphaPhi10);
    //     volScalarField alpha10("alpha10", alpha1);

    //     MULES::correct
    //     (
    //         geometricOneField(),
    //         alpha1,
    //         alphaPhiUn,
    //         talphaPhi1Corr.ref(),
    //         Sp,
    //         (-Sp*alpha1)(),
    //         oneField(),
    //         zeroField()
    //     );

    //     // Under-relax the correction for all but the 1st corrector
    //     if (aCorr == 0)
    //     {
    //         alphaPhi10 += talphaPhi1Corr();
    //     }
    //     else
    //     {
    //         alpha1 = 0.5*alpha1 + 0.5*alpha10;
    //         alphaPhi10 += 0.5*talphaPhi1Corr();
    //     }
    // }
    // else
    {

        volScalarField::Internal Su
        (
            fvc::ddt(rho1)*0.0
        );

        MULES::explicitSolve
        (
            geometricOneField(),
            alphaRho1,
            phi,
            alphaRhoPhi1,
            Su,
            Su,
            rho1.primitiveField(),
            zeroField()
        );
    }

    alpha1 = alphaRho1 / rho1;
    alpha2 = 1.0 - alpha1;

    mixture.correct();

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
