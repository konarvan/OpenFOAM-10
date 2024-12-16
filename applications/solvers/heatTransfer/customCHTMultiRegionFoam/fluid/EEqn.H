{
    volScalarField& he = thermo.he();

    volScalarField& T = thermo.T();
    volScalarField Cp = thermo.Cp();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, Cp * T)
      + (
            mvConvection.valid()
          ? mvConvection->fvmDiv(phi, Cp * T)
          : fvm::div(phi, Cp * T)
        )
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? mvConvection.valid()
            ? mvConvection->fvcDiv(fvc::absolute(phi, rho, U), p/rho)
            : fvc::div(fvc::absolute(phi, rho, U), p/rho)
          : -dpdt
        )
      + thermophysicalTransport.divq(Cp * T)
     ==
        rho*(U&g)
      + reaction.Qdot()
      + fvModels.source(rho, Cp * T)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    // Calculate entropy generation
    volScalarField Sgen
    (
        IOobject
        (
            "Sgen",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (thermophysicalTransport.divq(Cp * T) / T)
      + (turbulence.nuEff() / T) * magSqr(fvc::grad(U))
    );

    Sgen.write();

    fvConstraints.constrain(he);

    thermo.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}
