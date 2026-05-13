Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

IOdictionary FittingDict
(
    IOobject
    (
        "FittingDict",
        "constant",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

int numberOfMesh = FittingDict.lookupOrDefault<int>("numberOfMesh", 101);
scalar TMin        = FittingDict.lookupOrDefault<scalar>("minT", 300.0);
scalar TMax        = FittingDict.lookupOrDefault<scalar>("maxT", 3000.0);
scalar pressure    = FittingDict.lookupOrDefault<scalar>("pressure", 101325.0);

#include "createPreProcess.H"

PropertyReader propertyReader(mesh, speciesList);
propertyReader.generateTranCoeff(speciesList);

Info<< "Creating field binary diffusion coefficient Dij\n" << endl;
List<PtrList<volScalarField>> Dij(Y.size());

using MyThermoType = PPTransport<
    species::thermo<
        janafThermo<perfectGas<specie>>,
        sensibleEnthalpy
    >
>;

const auto& mixture =
    dynamic_cast<const DTMMultiComponentMixture<MyThermoType>&>
        (thermo.composition());

const PtrList<MyThermoType>& specieThermos = mixture.specieThermos();

forAll(TList, dataI)
{
    TLog[dataI] = Foam::log(TList[dataI]);
}

forAll(specieThermos, i)
{
    forAll(muList[i], celli)
    {
        scalar muVal     = specieThermos[i].mu(pList[celli], TList[celli]);
        scalar kappaVal  = specieThermos[i].kappa(pList[celli], TList[celli]);

        muList[i][celli]     = muVal;
        kappaList[i][celli]  = kappaVal;

        muLog[i][celli]      = Foam::log(muVal);
        kappaLog[i][celli]   = Foam::log(kappaVal);
    }
}

forAll(TList, celli)
{
    const auto& mixture_ = mixture.cellMixture(celli);

    forAll(DijList, i)
    {
        forAll(DijList[i], j)
        {
            scalar DijVal = mixture_.Dij(i, j, pList[celli], TList[celli]);
            DijList[i][j][celli] = DijVal;
            DijLog[i][j][celli]  = Foam::log(DijVal);
        }
    }

}

forAll(Y,i)
{
    Info << " =========== print speciesName ==============  " << endl;
    Info << " ith = " << i << ";--> speciesName =  " << composition.species()[i] << endl;

    // Fit Viscosity
    DPOLFT(N, TLog, muLog[i], W, MAXDEG, NDEG, EPS, R, IERR, muA);
    if (IERR == 1)
    {
        Info << "Successfully computed the polynomial for ." << endl;
    }
    else
    {
        Info << "Error occurred. IERR = " << IERR << endl;
    }
    DPCOEF(MAXDEG, C, muA, muTC[i]);

    Info << "Taylor series coefficients for mu:" << endl;

    forAll(muTC[i],index)
    {
        Info << "muTC(" << index << ") = " << muTC[i][index] << endl;
    }
    Info << endl;

    //Fit Thermal Conducvity
    DPOLFT(N, TLog, kappaLog[i], W, MAXDEG, NDEG, EPS, R, IERR, kappaA);
    if (IERR == 1)
    {
        Info << "Successfully computed the polynomial for ." << endl;
    }
    else
    {
        Info << "Error occurred. IERR = " << IERR << endl;
    }
    DPCOEF(MAXDEG, C, kappaA, kappaTC[i]);

    Info << "Taylor series coefficients for kappa:" << endl;

    forAll(kappaTC[i],index)
    {
        Info << "kappaTC(" << index << ") = " << kappaTC[i][index] << endl;
    }
    Info << endl;

    //Fit Binary Diffusion Coefficients
    forAll(speciesList, j)
    {
        DPOLFT(N, TLog, DijLog[i][j], W, MAXDEG, NDEG, EPS, R, IERR, DijA[j]);
        if (IERR == 1)
        {
            Info << "Successfully computed the polynomial for ." << endl;
        }
        else
        {
            Info << "Error occurred. IERR = " << IERR << endl;
        }
        DPCOEF(MAXDEG, C, DijA[j], DijTC[i][j]);

        Info << "Taylor series coefficients for Dij:" << endl;

        forAll(DijTC[i][j],index)
        {
            Info << "DijTC(" << index << ") = " << DijTC[i][j][index] << endl;
        }
        Info << endl;
    }

    Info << endl;

    propertyReader.printTranCoeff(runTime, composition.species()[i], muTC[i], kappaTC[i], DijTC[i]);


}
