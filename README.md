# Constant-Q-PSTD

This package includes four programs that can be used to do:

1) acoustic forward modeling

2) viscoacoustic forward modeling

3) acoustic reverse time migration (RTM)

4) Q-compensated reverse time migration (Q-RTM).

Both the forward modeling and RTM are based on pseudo-spectral time-domain (PSTD) numerical solutions of seismic wave equations. For the viscoacoustic wave modeling, the PSTD method is used to solve a fractional Laplacian constant-Q wave equation developed by Zhu and Harris (2014, Geophysics) and simplified by Chen et al. (2016, Geophysics). The optimal checkpointing technique (Symmes, 2007, Geophysics) is used in RTM and Q-RTM to avoid the storage of the source wavefields. 

All the codes are written with the standard C language. The MPI/GPU parallel programming techniques are adopted. Therefore, one should first install MPICH2 and CUDA C compilers correctly. Of course, before that, one should install a compatible NVIDIA driver to support the NVIDIA graphics cards.

#### Instructions to run the forward modeling packages ###

#1) Geometry generation
    a) Go to the folder "Geometry"
    b) Input: parameters.txt
    c) Compile and execute the codes by typing "Make" or "Make ge" in the command window
    d) Output: binary files in current path 
       tracenum.bin: total number of receivers for each shot
       shotxp.bin: x positions of all shots
       shotzp.bin: z positions (depth) of all shots
       gxp*.bin: the receiver's x position with * denoting the shot number 
(starting from 0)
	gzp*.bin: the receiver's z position with * denoting the shot number (starting from 0).

#2) Compiling the foward modeling codes
    a) Go to the folder "sgps-qmodeling"
    b) Specify the install path of the CUDA C compiler and MPICC compiler


To be continued






