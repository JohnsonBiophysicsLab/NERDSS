CXX ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -fopenmp -DNDEBUG
TARGET := reaction_propagation_omp_bench

.PHONY: all run clean
all: $(TARGET)

$(TARGET): reaction_propagation_omp_bench.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

run: $(TARGET)
	OMP_BENCH_THREADS=1,2,4,8 OMP_WAIT_POLICY=active OMP_PROC_BIND=spread OMP_PLACES=cores ./$(TARGET)

clean:
	rm -f $(TARGET)
