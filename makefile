CXX = g++
CXXFLAGS = -O3 -Wall -std=c++11
INCLUDE = -Iinclude
LDFLAGS = 
CUBA=/Users/jwalsh/Desktop/Research/HEP/Cuba-4.2/
GSLINC=/Users/jwalsh/Desktop/Research/HEP/gsl/include
GSLLIB=/Users/jwalsh/Desktop/Research/HEP/gsl/lib

# rule for object files
%.o: src/%.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) -I$(CUBA) -I$(GSLINC) $< -o $@ -L$(CUBA) -lcuba -L$(GSLLIB) -lgsl

# executables
all: test

test: test.o ThreeVector.o SPTkernels.o EFTkernels.o DiagramBase.o DiagramTree.o DiagramOneLoop.o Random.o Labels.o MomentumMap.o VertexMap.o Propagator.o LinearPowerSpectrumCAMB.o PowerSpectrum.o
	mkdir -p bin
	$(CXX) -o bin/$@ $^ $(CXXFLAGS) -L$(CUBA) -lcuba -L$(GSLLIB) -lgsl

clean:
	rm -f *.o
