#bin/sh
CUDA_INSTALL_PATH=/usr
#CUDA_INSTALL_PATH=/usr
LIBCUDA=-L$(CUDA_INSTALL_PATH)/lib64
LIB=-lcuda -lcudart -lcufft -lm
INC= -I $(CUDA_INSTALL_PATH)/include
#LIB2=/opt/apps/cuda/6.5/lib64/stubs/libcuda.so



# staggered-grid ps viscoacoustic modeling
vaps:visacousticsgps4.c func.o acsglw.o src.o
	mpicc -w -o vaps visacousticsgps4.c func.o acsglw.o src.o  $(LIBCUDA) $(LIB) $(INC) -lm -lstdc++ 
	rm -f *.o
	nohup mpirun -np 1 ./vaps &
func.o:func.c
	gcc -c func.c
sglowrank.o:sglowrank.c
	gcc -w -c sglowrank.c -lm
acsglw.o:acsglw.cu
	nvcc -c acsglw.cu -Wno-deprecated-gpu-targets
src.o : src.c 
	gcc -c -O3 src.c

# staggered-grid ps acoustic modeling
aps:acousticsgps4.c funca.o acps4.o src.o
	mpicc -w -o aps acousticsgps4.c funca.o acps4.o src.o  $(LIBCUDA) $(LIB) $(INC) -lm -lstdc++ -fPIE
	rm -f *.o
funca.o:funca.c
	gcc -c funca.c
acps4.o:acps4.cu
	nvcc -c acps4.cu -Wno-deprecated-gpu-targets
	nohup mpirun -np 1 ./aps &
	
	





