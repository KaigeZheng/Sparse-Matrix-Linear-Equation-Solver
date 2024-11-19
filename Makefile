# Compiler and Flags
CXX = g++
CXXFLAGS = -O2

# Lib Directory
INCLUDES = ./lib
DEMO = ./demo

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

# Clean
clean:
	rm -f test compute_det solver