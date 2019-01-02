#include "../code/Tri.hpp"

using namespace std;

int main(int argc, char* argv[]){


  ntuple<2, double>a({0, 0}), b({1, 0}), c({0, 1});

  Tri<2> tri(a, b, c);

  cout << "triangle :\n" << tri;


  return 0;
}
