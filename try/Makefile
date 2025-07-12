# File: Makefile

# Compiler and flags
NVCC        := nvcc
NVCC_FLAGS  := -O3 -arch=sm_70 -std=c++17

# Executable name
TARGET      := rk4_solver

# Source files
SRCS        := not.cu
OBJS        := $(SRCS:.cu=.o)

# Default rule
all: $(TARGET)

# Link rule
$(TARGET): $(OBJS)
	$(NVCC) $(NVCC_FLAGS) -o $@ $^

# Compile rule
%.o: %.cu
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)
