CXX = g++
CXXFLAGS = -std=c++11 -Wfatal-errors

MAIN_SRC = test_ntuple.cpp
TESTS = $(MAIN_SRC:.cpp=)
CODEREP = ../code/

all: $(TESTS)

test_ntuple: test_ntuple.o
test_ntuple.cpp: $(CODEREP)ntuple.hpp

test_fct: test_fct.o $(CODEREP)Function.o
	$(CXX) $(CXXFLAGS) $^ -o $@

test_fct_bad: test_fct_bad.o $(CODEREP)Function.o
	$(CXX) $(CXXFLAGS) $^ -o $@


%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@


clean:
	-rm *.o
	-rm *~


.PHONY: all clean
