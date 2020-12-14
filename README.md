# Constant-Q-PSTD

This package includes four programs that can be used to do:

1) acoustic forward modeling
2) viscoacoustic forward modeling
3) acoustic reverse time migration (RTM)
4) Q-compensated reverse time migration (Q-RTM).

Both the forward modeling and RTM are based on pseudo-spectral time-domain (PSTD) numerical solutions of seismic wave equations. For the viscoacoustic wave modeling, the PSTD method is used to solve a fractional Laplacian constant-Q wave equation developed by Zhu and Harris (2014, Geophysics) and simplified by Chen et al. (2016, Geophysics). The optimal checkpointing technique (Symmes, 2007, Geophysics) is used in RTM and Q-RTM to avoid the storage of the source wavefields. 

All the codes are written with the standard C language. The MPI/GPU parallel programming techniques are adopted. Therefore, one should first install MPICH2 and CUDA C compilers correctly. Of course, before that, one should install a compatible NVIDIA driver to support the NVIDIA graphics cards.

#### Instructions to run the packages ###

#1) Geometry generation

    a) Go to the folder "Geometry"
    b) Input: parameters.txt
    c) Compile and execute the codes by typing "Make" or "Make ge" in the command window
    d) Output: binary files in current path 
       tracenum.bin: total number of receivers for each shot
       shotxp.bin: x positions of all shots
       shotzp.bin: z positions (depth) of all shots
       gxp*.bin: the receiver's x position with * denoting the shot number (starting from 0)
       gzp*.bin: the receiver's z position with * denoting the shot number (starting from 0).

#2) Compile the foward modeling codes

    a) Go to the folder "sgps-qmodeling"
    b) Specify the install paths of the CUDA C compiler and MPICC compiler (if not installed in the default path) 
       CUDA_INSTALL_PATH=/usr
       LIBCUDA=-L$(CUDA_INSTALL_PATH)/lib64
       INC= -I $(CUDA_INSTALL_PATH)/include
       MPI_INSTALL_PATH=/usr/local
       LIBMPI=-L$(MPI_INSTALL_PATH)/lib
       INC2= -I $(MPI_INSTALL_PATH)/include
       LIB=-lcuda -lcudart -lcufft -lm
    c) Compiling commands: "Make aps" for acoustic PSTD codes and "Make vaps" for viscoacoustic PSTD codes
  
 #3) Execute the foward modeling codes
 
    a) Input parameters for acoustic forward modeling: go to the folder "input" and edit "ac2drealmodeling.txt" 
       especially fill in the following parameters correctly     
       --trace number of each shot in bin file (int)      // the output file in Setp #1)
       ../geometry/tracenum.bin
       --shot x position filename (float)                 // the output file in Setp #1)
       ../geometry/shotxp.bin
       --shot z position filenamee (float)                // the output file in Setp #1)
       ../geometry/shotzp.bin
       --receivers's x position prefix filename (float)   // the output file in Setp #1)
        ../geometry/gxp
       --receivers's z position prefix filename (float)   // the output file in Setp #1)
       ../geometry/gzp
       --input vp binary filename (float)         // velocity model stored in a binary format with x (trace id) as the primary column
       ../output/gasv160_398.bin
       --input density binary filename (float)    // density model stored in a binary format with x (trace id) as the primary column
       ../output/rho160_398.bin
       --output snapshots or not (1: yes, 0: no)  // if 1, ouput wavefield snapshot every 10 time steps to ../obsdata/snapa*_nz_nx.bin
         0
      b) For viscoacoustic forward modeling: go to the same folder "input" and edit "visac2drealmodeling.txt", similar to step a)
      c) Execute: "mpirun -np processnumber -machinefile hostlist ./aps" 
         Execute: "mpirun -np processnumber -machinefile hostlist ./vaps" 
         note that the processnumber should be specified and hostlist list the names of mpi nodes
         Or just typing "./aps" and "./vaps" to run the codes on the current node.
      d) Output: common-shot gathers in binary format
         "reca*_nr.bin" and "recva*_nr.bin" in the folder "obsdata"
         
   #4) Compile RTM and QRTM codes
   
     a) Go to the folder "sgps-qrtm"
     b) Specify the install paths of the CUDA C compiler and MPICC compiler (if not installed in the default path) 
       CUDA_INSTALL_PATH=/usr
       LIBCUDA=-L$(CUDA_INSTALL_PATH)/lib64
       INC= -I $(CUDA_INSTALL_PATH)/include
       MPI_INSTALL_PATH=/usr/local
       LIBMPI=-L$(MPI_INSTALL_PATH)/lib
       INC2= -I $(MPI_INSTALL_PATH)/include
       LIB=-lcuda -lcudart -lcufft -lm
     c) Compiling commands: "Make artm" for RTM and "Make qrtm" for QRTM
     
   #5) Execute RTM and QRTM codes
   
     a) Go to the folder "Input" to modify the input parameters: "a2drtm.txt" for RTM and "visac2drtm.txt" for QRTM
        especially notice the following parameters, taking "visac2drtm.txt" for example,
        
        --parameter to stabilize source ( x 1000000 )    // the actual value of sigma^2 is 1e-20 * 1e-6
          1e-20
        --parameter to stabilize receiver ( x 1000000 )  // the actual value sigma^2 is 1e-20 * 1e-6
          1e-20
        --trace number of each shot in bin file (int)    // the output file in Step #1)
          ../geometry/tracenum.bin
        --shot x position filename (float)               // the output file in Step #1)
          ../geometry/shotxp.bin
        --shot z position filenamee (float)               // the output file in Step #1)
          ../geometry/shotzp.bin
        --receivers's x position prefix filename (float)  // the output file in Step #1), but notice just giving the prefix name
          ../geometry/gxp
        --receivers's z position prefix filename (float)  // the output file in Step #1), but notice just giving the prefix name
          ../geometry/gzp
        --input vp binary filename (float)                // velocity model (binary format) for QRTM
          ../output/gasvm160_398.bin
        --input dens binary filename (float)              // density model (binary format) for QRTM
          ../output/rho160_398.bin
        --input Qp binary filename (float)                // Q model (binary format) for QRTM
          ../output/gasq160_398.bin
        --obs data binary data filename                   // observed data to be migrated, output file in Step #3), but notice just giving the prefix name
          ../obsdata/recva
          
        b) Execute: "mpirun -np processnumber -machinefile hostlist ./artm" for RTM 
         Execute: "mpirun -np processnumber -machinefile hostlist ./qrtm" for QRTM 
         note that the processnumber should be specified and hostlist list the names of mpi nodes
         Or just typing "./artm" and "./qrtm" to run the codes on the current node.
        c) Output: RTM profile stored in a binary file "RTM.BIN"
                   QRTM profile stored in a binary file "QRTM.BIN"
           Note that:  I)  The migration profiles have the same grid size (nx*nz) as the model defined in ../geometry/parameters.txt.
                       II) Additionally, each file (RTM.BIN or QRTM.BIN) includes two arrays with the same grid size of nx*nz.
                           The first one is the migration profile (denoted as r) and the second one is the source energy profile (s). 
                           Usually, we use the division of r by s (r./s) as the migration result.
                       
           
    
     
      

       




