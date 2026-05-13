Procedure:

#1. Convert Allclean and run.meshGeneration scripts into executable files: 
chmod +x run.meshGeneration 
chmod +x Allclean

#2. Clean the case file by running Allclean script:
./Allclean 

#3. Generate mesh using snappyHexMesh utility by running run.meshGeneration script:
./run.meshGeneration 

#4. Running simulations using catalystFoam solver:
catalystFoam 


# Note: 
   - You can run simulation using parallel mode since the mesh is quite large.
   - After generating mesh using snappyHexMesh, you can open the "system/controDict" 
   file and adjust the "writePrecision" from 14 to 6 if you want because such the 
   small precision is only neccessarly for snappyHexMesh. 

