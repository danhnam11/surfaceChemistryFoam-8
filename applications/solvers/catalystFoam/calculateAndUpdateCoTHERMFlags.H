label nSCoTHERM(CoTHERMMajorSpecies.size()); 

Info << "\nCalculating flags for coTHERM" << endl; 
// calculate flags for internal fields
forAll(flagSpecies, celli)
{
    //initialize the flag for species in a single cell 
    std::vector<double> flagSpeciesVector(nSCoTHERM, 1.0);  
    forAll(CoTHERMMajorSpecies, j)
    {

        if 
		(
			mag
			(
				Y[CoTHERMSpeciesIndex[j]].oldTime()[celli] 
			  - Y[CoTHERMSpeciesIndex[j]][celli]
			) <= epsilonSpecies
		)
        {
            flagSpeciesVector[j] = 0.0;
        }
        else
        {
            flagSpeciesVector[j] = 1.0;
        }
    }

    if 
	(
	    std::all_of
		(
		    flagSpeciesVector.begin(), 
			flagSpeciesVector.end(), 
			[&](double v){ return v == 0.0; }
		)
	)
    {
	    // Info << "flags of all species are equal " << endl;
	    flagSpecies[celli] = 0.0;
    }
    else
    {
	    // Info << "flags of all species are not equal: -> need to re-calculate Properties " << endl;	
	    flagSpecies[celli] = 1.0;
    }
}

// calculate flags for boundary fields
volScalarField::Boundary& pFlagSpecies = flagSpecies.boundaryFieldRef();

forAll(pFlagSpecies, patchi)
{
    forAll(pFlagSpecies[patchi], facei)
    {
        //initialize the flag for species in a single face 
        std::vector<double> flagSpeciesVectorF(nSCoTHERM, 1.0);
        forAll(CoTHERMMajorSpecies, j)
        {

            if 
            ( 
                mag
                (
                	Y[CoTHERMSpeciesIndex[j]].oldTime().boundaryField()[patchi][facei] 
                  - Y[CoTHERMSpeciesIndex[j]].boundaryField()[patchi][facei]
                ) <= epsilonSpecies
            )
            {
        	   flagSpeciesVectorF[j] = 0.0;
            }
            else
            {
        	   flagSpeciesVectorF[j] = 1.0;
            }
        }

        if 
		(
		    std::all_of
			(
			    flagSpeciesVectorF.begin(), 
				flagSpeciesVectorF.end(), 
				[&](double v){ return v == 0.0; }
			)
		)
        {
	        // Info << "flags of all species are equal " << endl;
	        pFlagSpecies[patchi][facei] = 0.0;
        }
        else
        {
	        // Info << "flags of all species are not equal: -> need to re-calculate Properties " << endl;	
	        pFlagSpecies[patchi][facei] = 1.0;
        }
       
    }
}

// update flags 
thermo.updateCoTHERMFlags(flagSpecies); 
Info << "End of calculating flags for coTHERM\n" << endl; 
