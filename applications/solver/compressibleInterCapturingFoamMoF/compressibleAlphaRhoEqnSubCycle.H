
#include "rhofs.H"
// tmp<surfaceScalarField> talphaPhi1(alphaPhi10);

// if (nAlphaSubCycles > 1)
// {
//     dimensionedScalar totalDeltaT = runTime.deltaT();

//     talphaPhi1 = new surfaceScalarField
//     (
//         IOobject
//         (
//             "alphaPhi1",
//             runTime.timeName(),
//             mesh
//         ),
//         mesh,
//         dimensionedScalar(alphaPhi10.dimensions(), Zero)
//     );

//     surfaceScalarField rhoPhiSum
//     (
//         IOobject
//         (
//             "rhoPhiSum",
//             runTime.timeName(),
//             mesh
//         ),
//         mesh,
//         dimensionedScalar(rhoPhi.dimensions(), Zero)
//     );

//     tmp<volScalarField> trSubDeltaT;

//     if (LTS)
//     {
//         trSubDeltaT =
//             fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
//     }

//     for
//     (
//         subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
//         !(++alphaSubCycle).end();
//     )
//     {
//         #include "alphaRhoEqn.H"
//         talphaPhi1.ref() += (runTime.deltaT()/totalDeltaT)*alphaPhi10;
//         rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
//     }

//     rhoPhi = rhoPhiSum;
// }
// else
{
    #include "alphaRhoEqn.H"
}

alphaRho2 = alpha2*rho2;


alphaPhi1 = alphaRhoPhi1 / rho1f;
alphaPhi2 =  phi - alphaPhi1;
alphaRhoPhi2 = (phi - alphaPhi1)*rho2f;


rhoPhi = alphaRhoPhi1 + alphaRhoPhi2;
rho == alphaRho1 + alphaRho2;

surfaceScalarField alphaf = fvc::interpolate(alpha1);
forAll(alphaf,facei)
{
    if (phi[facei] != 0)
    {
        alphaf[facei] = alphaPhi1[facei]/phi[facei];
    }
}

rhof = alphaf*rho1f + (1-alphaf)*rho2f;


volScalarField::Internal contErr
(
    (
        fvc::ddt(rho) + fvc::div(rhoPhi)
      - (fvOptions(alpha1, mixture.thermo1().rho())&rho1)
      - (fvOptions(alpha2, mixture.thermo2().rho())&rho2)
    )()
);
