#bin/sh
CUDA_INSTALL_PATH=/home/chenhm/programs/cuda
#CUDA_INSTALL_PATH=/usr
LIBCUDA=-L$(CUDA_INSTALL_PATH)/lib64
LIB=-lcuda -lcudart -lcufft -lm
INC= -I $(CUDA_INSTALL_PATH)/include
#LIB2=/opt/apps/cuda/6.5/lib64/stubs/libcuda.so



# staggered-grid visacoustic wave rtm
vaps:visacousticsgps4.c func.o acps4.o src.o
	mpicc -w -o vaps visacousticsgps4.c func.o acps4.o src.o  $(LIBCUDA) $(LIB) $(INC) -lm
func.o:func.c
	gcc -c func.c
acps4.o:acps4.cu
	nvcc -c acps4.cu -Wno-deprecated-gpu-targets
src.o : src.c 
	gcc -c -O3 src.c




