CXX      ?= g++
CXXFLAGS ?= -O2 -fopenmp
SRC      := $(wildcard src/*.cpp)

secret: $(SRC) $(wildcard src/*.h)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC) -lz

clean:
	rm -f secret
