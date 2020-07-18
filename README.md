# code-for-paper-fracture-solute-transport
This code is for solute transport in  rough fractures!
The file named "geoSMOOTH.dat" is one of the example geometry data.
Code running steps:            
(1) Based on Linux system, download and unzip the palabos package. Download address: https://github.com/gladk/palabos            
(2) Copy the folder "code for paper" into the "example" folder in the directory of palabos,and create a new empty file named "tmp" to restore the simulation results.            
(3) Open the terminal and enter "make";            
(4) After make is completed, enter "mpirun - np 32 ./ SoluteTransport", where 32 threads are used for parallel operation.
this simulation may cost long time which is more than 2 days.
