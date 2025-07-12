#path: /usr/local/cuda/bin
#    nvcc -L /usr/local/sundials/lib64 
all:
	nvcc --resource-usage startup.cu -o startup.x -Xcompiler -fopenmp

clean:
	rm -f *~ *.x
