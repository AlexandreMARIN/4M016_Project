CXX = g++
CXXFLAGS = -std=c++11 -Wfatal-errors

SRC = Function.cpp Tri.cpp ItgQuadForm.cpp
HEADERS = ntuple.hpp Matrix.hpp

OBJ = $(SRC:.cpp=.o)

all: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o project



Tri.o: ntuple.hpp

ItgQuadForm.o: ntuple.hpp Matrix.hpp ItgQuadForm.cpp ItgQuadForm.hpp
	$(CXX) $(CXXFLAGS) -DDIR_COORDS=\"\" -c $< -o $@

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@



clean:
	-rm *.o
	-rm *~


.PHONY: all clean
