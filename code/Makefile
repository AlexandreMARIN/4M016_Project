CXX = g++
CXXFLAGS = -std=c++11 -Wfatal-errors

SRC = Function.cpp

OBJ = $(SRC:.cpp=.o)

all: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o project


Function.cpp: Function.hpp

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@



clean:
	-rm *.o
	-rm *~


.PHONY: all clean
