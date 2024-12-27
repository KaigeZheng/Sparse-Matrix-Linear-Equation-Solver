# Compiler and Flags
CXX = g++
CXXFLAGS = -O2

# Lib Directory
INCLUDES = ./lib
DEMO = ./demo
EIGEN = ~/eigen-3.4.0/

# Default Target
all: test

# test Target
test: test.cpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDES) test.cpp -o test

# compute_det Target
compute_det: $(DEMO)/compute_det.cpp
	$(CXX) $(CXXFLAGS) $(DEMO)/compute_det.cpp -o compute_det

# solver Target
solver: $(DEMO)/solver.cpp
	$(CXX) $(CXXFLAGS) $(DEMO)/solver.cpp -o solver

baseline:baseline.cpp
	$(CXX) baseline.cpp -o baseline

eigen:eigen-optimized.cpp
	$(CXX) eigen-optimized.cpp -I$(EIGEN)-o eigen-optimized

best:CSR-optimized.cpp
	$(CXX) CSR-optimized.cpp -o CSR-best

simd:CSR-simd.cpp
	$(CXX) CSR-simd.cpp -o CSR-simd

# Clean
clean:
	rm -f test compute_det solver baseline eigen-optimized CSR-best CSR-simd