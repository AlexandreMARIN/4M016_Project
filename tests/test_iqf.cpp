#include "../code/ItgQuadForm.hpp"

using namespace std;
using namespace fct;

int main(int argc, char* argv[]){

  Tri<2> tri;
  tri[0][0] = 0.;
  tri[0][1] = 0.;
  tri[1][0] = 1.;
  tri[1][1] = 0.;
  tri[2][0] = 0.;
  tri[2][1] = 1.;

  cout << "linear system from the right triangle (don't take into account the last column):\n";
  cout << ItgQuadForm::get_lin_sys(tri);

  cout << "\npoints/weights for exact integration of polynomials of degree 8:\n";
  const vector<ntuple<3, double> >& pw8 = ItgQuadForm::getpw(8);

  for(int i=0;i<pw8.size();i++){
    cout << pw8[i] << "\n";
  }

  return 0;
}
