CXX = g++
CXXFLAGS = -std=c++11 -Wfatal-errors

MAIN_SRC = test_ntuple.cpp test_fct.cpp test_fct_bad.cpp test_tri.cpp test_mat.cpp test_iqf.cpp test_integ.cpp
TESTS = $(MAIN_SRC:.cpp=)
DIRCODE = ../code/

all: $(TESTS)

test_ntuple: test_ntuple.o $(DIRCODE)ntuple.hpp
	$(CXX) $(CXXFLAGS) test_ntuple.o -o $@

test_fct: test_fct.o $(DIRCODE)Function.o $(DIRCODE)Tri.o $(DIRCODE)ItgQuadForm.o
	$(CXX) $(CXXFLAGS) $^ -o $@

test_fct_bad: test_fct_bad.o $(DIRCODE)Function.o $(DIRCODE)Tri.o $(DIRCODE)ItgQuadForm.o
	$(CXX) $(CXXFLAGS) $^ -o $@

test_tri: test_tri.o $(DIRCODE)Tri.o $(DIRCODE)Tri.hpp
	$(CXX) $(CXXFLAGS) test_tri.o $(DIRCODE)Tri.o -o $@

test_mat: test_mat.o $(DIRCODE)Matrix.hpp
	$(CXX) $(CXXFLAGS) test_mat.o -o $@

test_iqf: test_iqf.o $(DIRCODE)ItgQuadForm.o $(DIRCODE)Tri.o $(DIRCODE)Tri.hpp
	$(CXX) $(CXXFLAGS) $^ -o $@

$(DIRCODE)ItgQuadForm.o: $(DIRCODE)ItgQuadForm.cpp $(DIRCODE)Tri.hpp
	$(CXX) $(CXXFLAGS) -DDIR_COORDS=\"../code/\" -c $< -o $@

test_integ: test_integ.o $(DIRCODE)ItgQuadForm.o $(DIRCODE)Tri.o $(DIRCODE)Function.o $(DIRCODE)Matrix.hpp $(DIRCODE)Tri.hpp
	$(CXX) $(CXXFLAGS) $^ -o $@


%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	-rm *.o
	-rm *~


.PHONY: all clean
