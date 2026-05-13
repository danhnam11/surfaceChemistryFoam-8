Info<< "Reading thermophysical properties for catalyst\n" << endl;
autoPtr<psiReactionCatThermo> pCatThermo(psiReactionCatThermo::New(mesh));
psiReactionCatThermo& catThermo = pCatThermo();
catThermo.validate(args.executable(), "h", "e");
Info<< "End of creating thermophysical properties for catalyst\n" << endl;

Info<< "\nCreating surface chemsitry model for catalyst\n" << endl;
autoPtr<BasicSurfaceChemistryModel<psiReactionThermo, psiReactionCatThermo>> 
surfaceReaction
(   
    BasicSurfaceChemistryModel<psiReactionThermo, psiReactionCatThermo>::
    New(thermo, catThermo)
);
Info<< "\nEnd of creating surface chemistry model for catalyst\n" << endl;
