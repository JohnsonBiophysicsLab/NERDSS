NVCC ?= nvcc
COMPUTE_CAP ?= $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1)
GPU_ARCH ?= sm_$(subst .,,$(COMPUTE_CAP))
NVCCFLAGS ?= -O3 -std=c++17 -arch=$(GPU_ARCH)
TARGET := reaction_propagation_bench

.PHONY: all run clean
all: $(TARGET)

$(TARGET): reaction_propagation_bench.cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)
