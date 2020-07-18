# code-for-paper-fracture-solute-transport
this code is for solute transport in  rough fractures!

Code running steps:            
(1) Based on Linux system, download and unzip the palabos package. Download address: https://github.com/gladk/palabos            
(2) Copy the folder "code for paper" into the "example" folder in the directory of palabos.            
(3) Open the terminal and enter "make";            
(4) After make is completed, enter "mpirun - np 32 ./ SoluteTransport", where 32 threads are used for parallel operation.
this simulation may cost long time which is more than 2 days.
