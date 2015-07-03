CXX = g++
CXXFLAGS = -O3 -Wall -std=c++11
INCLUDE = -Iinclude
LDFLAGS = 
CUBA=/Users/jwalsh/Desktop/Research/HEP/Cuba-4.2/

# rule for object files
%.o: src/%.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -I$(CUBA) $< -o $@ -L$(CUBA) -lcuba

# executables
# KernelBase.o 
test: test.o SPTkernels.o EFTkernels.o Diagram.o Random.o Labels.o Momentum.o Propagator.o PowerSpectrum.o Bispectrum.o Trispectrum.o
	mkdir -p bin
	$(CXX) -o bin/$@ $^ $(CXXFLAGS) -L$(CUBA) -lcuba

clean:
	rm -f *.o
