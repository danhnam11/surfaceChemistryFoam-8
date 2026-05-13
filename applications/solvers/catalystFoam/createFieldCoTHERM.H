
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

bool usingDetailedTransportModel
(
    thermDict.lookupOrDefault("usingDetailedTransportModel", false)
); 

bool usingCoTHERM
(
    thermDict.lookupOrDefault("usingCoTHERM", false)
); 

bool usingCoTHERMOnlyT
(
    thermDict.lookupOrDefault("usingCoTHERMOnlyT", false)
); 

bool soretEffect
(
    (thermDict.found("soretEffects"))
    ? thermDict.lookupOrDefault("soretEffects", false)
    :thermDict.lookupOrDefault("soretEffect", false)
); 

if (soretEffect)
{
    Info << "\nSoret effect is considered" << endl;
}

//List<word> CoTHERMMajorSpecies(thermDict.lookup("majorSpeciesForCoTHERM"));
List<word> CoTHERMMajorSpecies
(
    (usingCoTHERM || usingCoTHERMOnlyT)
  ? List<word>(thermDict.lookup("majorSpeciesForCoTHERM"))
  : List<word>(1, inertSpecie)
);
scalar epsilonT(thermDict.lookupOrDefault("epsilonT", 0.1));
scalar epsilonSpecies(thermDict.lookupOrDefault("epsilonS", 0.001));

if (usingCoTHERM || usingCoTHERMOnlyT)
{
    Info << "\nMajor species are selected for CoTHERM are " 
         << CoTHERMMajorSpecies << endl;
}
else
{
    Info << "\ncoTHERM is not used" << endl; 
    Info << "inertSpecies is temporarily selected" 
         << "for coTHERM (eventhough it is not used)" << endl;  
    Info << "\nMajor species are selected for CoTHERM are " 
         << CoTHERMMajorSpecies << endl;    
}
Info << "epsilon for Temperature in CoTHERM = " << epsilonT << " [K]" << endl;
Info << "epsilon for Species in CoTHERM     = " 
     << epsilonSpecies << " [-](in mass fraction)\n" << endl;

List<label> CoTHERMSpeciesIndex(CoTHERMMajorSpecies.size());

forAll(CoTHERMSpeciesIndex, i)
{
    if (!composition.species().found(CoTHERMMajorSpecies[i]))
    {
        Info << "\nError!, cannot find " << CoTHERMMajorSpecies[i] 
             << " in the chemical Mechanism \n" << endl;
        Info << "please adjust the major species list" 
             << " in thermophysicalProperties dictionary!"  
             << exit(FatalIOError);
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
