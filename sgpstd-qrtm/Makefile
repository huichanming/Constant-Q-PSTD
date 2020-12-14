#bin/sh
CUDA_INSTALL_PATH=/usr
#CUDA_INSTALL_PATH=/usr
LIBCUDA=-L$(CUDA_INSTALL_PATH)/lib64
LIB=-lcuda -lcudart -lcufft -lm
INC= -I $(CUDA_INSTALL_PATH)/include
#LIB2=/opt/apps/cuda/6.5/lib64/stubs/libcuda.so



# staggered-grid lowrank acoustic wave rtm
qrtm:acrtmain.c func.o acsglw.o src.o
	mpicc -w -o qrtm acrtmain.c func.o acsglw.o src.o  $(LIBCUDA) $(LIB) $(INC) -lm -lstdc++ 
	rm -f *.o
	nohup mpirun -np 1 ./qrtm &
func.o:func.c
	gcc -c func.c
sglowrank.o:sglowrank.c
	gcc -w -c sglowrank.c -lm
acsglw.o:acsglw.cu
	nvcc -c acsglw.cu -Wno-deprecated-gpu-targets
src.o : src.c 
	gcc -c -O3 src.c


artm:main.c funca.o acps4.o src.o
	mpicc -w -o artm main.c funca.o acps4.o src.o  $(LIBCUDA) $(LIB) $(INC) -lm -lstdc++ 
	rm -f *.o
	nohup mpirun -np 1 ./artm &
funca.o:funca.c
	gcc -c funca.c
acps4.o:acps4.cu
	nvcc -c acps4.cu -Wno-deprecated-gpu-targets
	
	





