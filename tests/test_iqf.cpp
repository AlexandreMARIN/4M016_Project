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


  vector<ntuple<3, double> > pw;
  for(int d : vector<int>({2, 4, 5, 7, 9, 11, 13, 14, 16, 18, 20, 21, 23, 25})){
    pw = ItgQuadForm::getpw(d);
    cout << "\npoints/weights for exact integration of polynomials of degree "<< d <<":\n";
    cout << "number N of points: " << pw.size() << "\n";
    for(int i=0;i<pw.size();i++){
      cout << pw[i] << "\n";
    }
  }


  return 0;
}
