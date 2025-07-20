SUNDIALS_DIR := /usr/local/sundials
INCLUDE_FLAGS := -I$(SUNDIALS_DIR)/include
LIBRARY_FLAGS := -L$(SUNDIALS_DIR)/lib64

LDLIBS := -lsundials_cvode -lsundials_nveccuda \
          -lsundials_sunlinsolcusolversp -lsundials_sunmatrixcusparse \
          -lsundials_core -lcusolver -lcusparse -lcudart -lcudadevrt -lm

NVCC := nvcc
CUDA_ARCH := -arch=sm_89
CUDA_DEBUG := -G

TARGET := 3dODE_gpu
SRC := 3d_ode.cu

all: $(TARGET)

$(TARGET): $(SRC)
	$(NVCC) $(CUDA_ARCH) $(CUDA_DEBUG) $(INCLUDE_FLAGS) -o $@ $< $(LIBRARY_FLAGS) $(LDLIBS)

run: $(TARGET)
	LD_LIBRARY_PATH=$(SUNDIALS_DIR)/lib64 ./$(TARGET)

clean:
	rm -f $(TARGET)
