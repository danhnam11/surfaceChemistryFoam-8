# surfaceChemistryFoam-8

## General Information
A new surface reaction library that extends OpenFOAM(OF) with robust computational fluid dynamics (CFD) capabilities for simulating reactive flows coupled with detailed surface chemical kinetics. The library supports multiple reaction formulations, including Arrhenius-type expressions, sticking coefficient models, and surface coverage-dependent mechanisms. The developed library is native OF code and can be easily integrated into any OF-based solver (of the same version) to simulate catalytic reacting flows, supporting both particle-resolved approache. Readers are referred to our paper for all validation data.

## Installation
- The complete installation of the OpenFOAM 8.0 framework in a Linux operating system is required before installing this package, as it is designed for the Linux-based OpenFOAM 8.0 version. 
- Prepare a directory on your system, for example, _yourDirectory_:

		mkdir ~/OpenFOAM/yourDirectory/
		cd ~/OpenFOAM/yourDirectory/	
- Download source files using git: 

		git clone https://github.com/danhnam11/surfaceChemistryFoam-8.git

- Specify the path of the _src_ directory of this package to an environment variable named _LIB_SCHEM8_SRC_. Suppose the _surfaceChemistryFoam-8_ have been downloaded into _yourDirectory_. Then the following commands should be executed to specify the path of the _src_:

		echo "export LIB_SCHEM8_SRC=~/OpenFOAM/yourDirectory/surfaceChemistryFoam-8/src/" >> ~/.bashrc
		source ~/.bashrc

- To compile the necessary libraries and solver, go to _surfaceChemistryFoam-8_ directory and run the _Allwmake_ script (it may take hours to be finished):

		cd ~/OpenFOAM/yourDirectory/surfaceChemistryFoam-8/
		./Allwmake

- After successful compilation, the following libraries are saved at _$FOAM_USER_LIBBIN_ :

		libcatSpecie.so
		libthermophysicalProperties.so
		libfluidThermophysicalModels.so	
		libreactionThermophysicalModels.so
		libcatReactionThermophysicalModels.so		
		libchemistryModel.so		
		libsurfaceChemistryModel.so		
		libsolidThermo.so
		libSLGThermo.so
		libfluidThermoMomentumTransportModels.so
		libthermophysicalTransportModels.so
		libpsiReactionThermophysicalTransportModels.so
		librhoReactionThermophysicalTransportModels.so
		libradiationModels.so
		libcombustionModels.so
		libfvOptions.so
		libspecieTransfer.so		
		libcatalyticSurfaceFvPatchFields.so
		libsurfaceFilmModels.so
		libsurfaceFilmDerivedFvPatchFields.so
		libthermalBaffleModels.so
		libfieldFunctionObjects.so
		libforces.so
		liblagrangianFunctionObjects.so
		libsolverFunctionObjects.so
		libutilityFunctionObjects.so
		
- and the following executable programs are saved at _$FOAM_USER_APPBIN_ :

		catalystFoam
		surfChemkinToFoam
		
		DTLreactingFoam
		DTMchemkinToFoam
		FTMchemkinToFoam
		
		reactingFoam			
		chemkinToFoam

- These newly compiled libraries, solvers, and utilities are now ready for use.
- It is important to note that if a different solver (i.e., program) relies on any of the aforementioned compiled libraries, its corresponding _options_ file, located in the _Make_ directory, must be updated accordingly. The solver should then be recompiled to prevent potential conflicts, such as segmentation faults. For reference, the _Make_ directory of the _reactingFoam_ solver included in this package provides a convenient example. Although this version of _reactingFoam_ is nearly identical to the original, its _options_ file has been modified to link against compiled libraries in this package. As a result, it can be used seamlessly in place of the original version without issue.

- To remove all compiled libraries and solvers, go to _surfaceChemistryFoam-8_ directory and run the _Allwclean_ script:

		cd ~/OpenFOAM/yourDirectory/surfaceChemistryFoam-8/
		./Allwclean

## Using this package 
Upon completing the compilation process, the _catalystFoam_ solver can be utilized by simply typing its name in the terminal, _catalystFoam_. All important instructions for using surface chemistry library in a new developed OF-based solver and case setting are provided in _documentations_ directory.

## Tutorials
A test cases in catalytic processes is available in the _tutorials_ directory.

	cd ~/OpenFOAM/yourDirectory/surfaceChemistryFoam-8/tutorials/

## Authors 
This package was developed based on the _DTLreactingFoam_ package [1] at the Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering, Ulsan National Institute of Science and Technology (UNIST), Korea (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/). If you publish results obtained by using this package, please cite our paper as follows:
- D. N. Nguyen, J. H. Lee, H. W. Seo, H. J. Ahn, C. S. Yoo, surfaceChemistryFoam: An OpenFOAM-based library for detailed surface chemistry in reacting flow simulations, Computer Physics Communications (2026)(submitted).

- D. N. Nguyen, J. H. Lee, C. S. Yoo, DTLreactingFoam: An efficient CFD tool for laminar reacting flow simulations using detailed chemistry and transport with time-correlated thermophysical properties, Computer Physics Communications 322 (2026) 110052 (https://doi.org/10.1016/j.cpc.2026.110052).

Contact:
- danhnam11@gmail.com or nam.nguyendanh@hust.edu.vn 

## Reference
[1] D. N. Nguyen, J. H. Lee, C. S. Yoo, DTLreactingFoam: An efficient CFD tool for laminar reacting flow simulations using detailed chemistry and transport with time-correlated thermophysical properties, Computer Physics Communications 322 (2026) 110052 (https://doi.org/10.1016/j.cpc.2026.110052).
