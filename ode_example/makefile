rm CUDAFunction.o
rm CUDAInterface.o
rm kernel.o
rm CHAMR.o
rm HAMR.o
rm point.o
rm polygon_io.o
rm reciprocity.o
rm Summary.o
rm waveform.o
rm HAMR

/usr/local/cuda/bin/nvcc -w -m64 -O2 -I/usr/local/cuda/include -std=c++17 -c HAMR.cpp
/usr/local/cuda/bin/nvcc -w -m64 -O2 -I/usr/local/cuda/include -std=c++17 -c CHAMR.cpp

/usr/local/cuda/bin/nvcc -w -m64 -O2 -I/usr/local/cuda/include -std=c++17 -c point.cpp
/usr/local/cuda/bin/nvcc -w -m64 -O2 -I/usr/local/cuda/include -std=c++17 -c polygon_io.cpp
/usr/local/cuda/bin/nvcc -w -m64 -O2 -I/usr/local/cuda/include -std=c++17 -c reciprocity.cpp
/usr/local/cuda/bin/nvcc -w -m64 -O2 -I/usr/local/cuda/include -std=c++17 -c waveform.cpp
/usr/local/cuda/bin/nvcc -w -m64 -O2 -I/usr/local/cuda/include -std=c++17 -c Summary.cpp

/usr/local/cuda/bin/nvcc -w -m64 -O2 -arch=compute_75 -code=compute_75,sm_75 -c CUDAFunction.cu
/usr/local/cuda/bin/nvcc -w -m64 -O2 -arch=compute_75 -code=compute_75,sm_75 -c CUDAInterface.cu
/usr/local/cuda/bin/nvcc -w -m64 -O2 -arch=compute_75 -code=compute_75,sm_75 -c kernel.cu

/usr/local/cuda/bin/nvcc -w -m64 -O2 -arch=compute_75 -code=compute_75,sm_75 -std=c++17 -lcudart -L/usr/local/cuda/lib64  CUDAFunction.o CUDAInterface.o kernel.o CHAMR.o HAMR.o point.o polygon_io.o reciprocity.o Summary.o waveform.o -o HAMR


