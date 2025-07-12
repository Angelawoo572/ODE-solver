APPS=hello

all: ${APPS}

%: %.cu
	nvcc -O2 -arch=sm_89 -o $@ $<
clean:
	rm -f ${APPS}
