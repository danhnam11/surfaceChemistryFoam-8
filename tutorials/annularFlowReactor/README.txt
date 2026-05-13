#step 1. generate mesh using blockMesh by running
blockMesh 

#step 2. running simulation by typing catalystFoam in the ternimal 
catalystFoam

#step 3. post-processing to get O2 conversion 
#After running simulation, "postProcessing" directory will be generated.
#To get O2 conversion, run the postO2Conversion executable file by typing in the ternimal
./postO2Conversion 

# in case the postO2Conversion is not executable file, you can change postO2Conversion into 
# the executable file by using 
chmod +X postO2Conversion 

# or you can compile the postO2Conversion.c file, then run it 

# After running postO2Conversion, there is something on the screen like:
time_O2       = 0.000286394
time_N2       = 0.000286394
XO2out        = 9.4340700000e-03
XN2out        = 9.5055320000e-01
O2Conversion  = 5.714204 %

# so the last line is the O2 conversion result.

