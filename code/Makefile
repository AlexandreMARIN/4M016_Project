CXX = g++
CXXFLAGS = -std=c++11 -Wfatal-errors

SRC = Function.cpp Tri.cpp ItgQuadForm.cpp project.cpp
TEMPL = ntuple.hpp Matrix.hpp Tri.hpp

OBJ = $(SRC:.cpp=.o)

all: $(OBJ) $(TEMPL)
	$(CXX) $(CXXFLAGS) $(OBJ) -o project

project.o: $(TEMPL) ItgQuadForm.hpp Function.hpp

Tri.o: ntuple.hpp Tri.hpp

ItgQuadForm.o: $(TEMPL) ItgQuadForm.hpp
	$(CXX) $(CXXFLAGS) -DDIR_COORDS=\"\" -c ItgQuadForm.cpp -o $@

Function.o: Function.hpp $(TEMPL) ItgQuadForm.hpp

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $*.cpp -o $@



clean:
	-rm *.o
	-rm *~


.PHONY: all clean
