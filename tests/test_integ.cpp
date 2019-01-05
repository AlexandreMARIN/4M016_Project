#include <iostream>

#include "../code/Function.hpp"

using namespace std;
using namespace fct;

double one(double x, double y){
  return 1.0;
}

int main(int argc, char* argv[]){

  R2toRCppFct One(&one);
  LinCbnFct myfct({&pr2, &One}, {5.4, 0.6});
  ProdFct xy({&pr1, &pr2});
  ntuple<2, double>a({0, 0}), b({1, 0}), c({0, 1});
  Tri<2> tri(a, b, c);
  ntuple<3, double> proj;

  cout << "K is the triangle: " << tri;
  cout << "\nintegral of pr1 on K :\n" << pr1.integrate(tri);

  cout << "\nintegral on K of (x, y)|->xy:\n" << xy.integrate(tri) << "\n\n";

  cout << "projection on P1(K) of pr1:\n" << pr1.projection_on_P1(tri);
  cout << "\n\nprojection on P1(K) of myfct:\n" << myfct.projection_on_P1(tri) << "\n\n";

  cout << "E_K for myfct:\n" << myfct.get_E_K(tri, proj) << "\n\n";

  setepsilon(1e-7);
  xy.exportGnuplot(Tri<2>::generateMesh(-1., 1., -1., 1., 4), "out_test.gp");

  return 0;
}
