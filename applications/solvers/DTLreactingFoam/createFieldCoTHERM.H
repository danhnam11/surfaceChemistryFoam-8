
IOdictionary thermDict
(
    IOobject
    (
        "thermophysicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

List<word> CoTHERMMajorSpecies(thermDict.lookup("majorSpeciesForCoTHERM"));
scalar epsilonT(thermDict.lookupOrDefault("epsilonT", 0.1));
scalar epsilonSpecies(thermDict.lookupOrDefault("epsilonS", 0.001));

Info << "\nMajor species are selected for CoTHERM are " << CoTHERMMajorSpecies << endl;
Info << "epsilon for Temperature in CoTHERM = " << epsilonT << " [K]" << endl;
Info << "epsilon for Species in CoTHERM     = " << epsilonSpecies << " [-](in mass fraction)\n" << endl;

List<label> CoTHERMSpeciesIndex(CoTHERMMajorSpecies.size());

forAll(CoTHERMSpeciesIndex, i)
{
    if (!composition.species().found(CoTHERMMajorSpecies[i]))
    {
        Info << "\nError!, cannot find " << CoTHERMMajorSpecies[i] << " in the chemical Mechanism \n" << endl;
        Info << "please adjust the major species list in coTranProperites dictionary!"  << exit(FatalIOError);
    }
    else
    {
        CoTHERMSpeciesIndex[i] = composition.species()[CoTHERMMajorSpecies[i]];
    }
}

volScalarField flagSpecies
(
    IOobject
    (
        "flagSpecies",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("flagSpecies", dimensionSet(0, 0, 0, 0, 0), 1.0)
);
