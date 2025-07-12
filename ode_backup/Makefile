# CUDA 编译器
NVCC = nvcc

# 编译选项（为 sm_89 优化并 suppress 警告）
NVCC_FLAGS = -O3 -gencode arch=compute_89,code=sm_89 -Wno-deprecated-gpu-targets

# 注意路径：从 ode_backup/ 目录出发
SRCS = main.cu ../culsoda/trunk/Code/cuLsoda.cu ../culsoda/trunk/Code/cuLsoda_kernel.cu

# 可执行程序名
TARGET = culsoda

# 默认目标
all: $(TARGET)

$(TARGET): $(SRCS)
	$(NVCC) $(NVCC_FLAGS) $(SRCS) -o $(TARGET)

clean:
	rm -f $(TARGET) *.o *.out
