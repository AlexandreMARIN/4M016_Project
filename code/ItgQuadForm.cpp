#include "ItgQuadForm.hpp"

#include <cmath>

//#define DIR_COORDS "./"

using namespace std;
using namespace fct;

Matrix<3, 4> ItgQuadForm::mat;
map<int, vector<ntuple<3, double> > > ItgQuadForm::pts_w;
const ItgQuadForm* const ItgQuadForm::obj = new ItgQuadForm;

ItgQuadForm::ItgQuadForm(){

  string filename{DIR_COORDS};
  char buf[100];
  int deg, nb, p;
  vector<ntuple<3, double> > pw;

  filename += "coords.txt";

  ifstream coords{filename};

  while(coords){
    coords.read(buf, 19);//integration degree=
    coords >> deg;
    coords.read(buf, 4);//__N=
    coords >> nb;
    coords.read(buf, 1);//:
    pw.resize(nb);
    for(p=0;p<nb;p++){
      coords >> pw[p][0] >> pw[p][1] >> pw[p][2];//two coordinates, one weight
    }
    pts_w[deg] = pw;
    coords.read(buf, 1);
  }

  coords.close();

}

const std::vector<ntuple<3, double> >& ItgQuadForm::getpw(int d){
  auto iter = pts_w.lower_bound(d);

  if(iter==pts_w.end()){
    throw out_of_range("ItgQuadForm[] : bad argument\n");
  }

  return iter->second;
}

Matrix<3, 4>& ItgQuadForm::get_lin_sys(const Tri<2>& tri){

  double area = 0.5*abs((tri[1][0]-tri[0][0])*(tri[2][1]-tri[0][1]) - (tri[1][1]-tri[0][1])*(tri[2][0]-tri[0][0]));


  mat(0, 2) = area;

  mat(0, 0) = area*(tri[0][0]+tri[1][0]+tri[2][0])/3.0;
  mat(1, 2) = mat(0, 0);

  mat(0, 1) = area*(tri[0][1]+tri[1][1]+tri[2][1])/3.0;
  mat(2, 2) = mat(0, 1);

  mat(1, 1) = area*(tri[0][0]*(2.*tri[0][1]+tri[1][1]+tri[2][1])
		   +tri[1][0]*(tri[0][1]+2.*tri[1][1]+tri[2][1])
		   +tri[2][0]*(tri[0][1]+tri[1][1]+2.*tri[2][1]))/12.;
  mat(2, 0) = mat(1, 1);

  mat(1, 0) = area*(tri[0][0]*tri[0][0] + tri[0][0]*tri[1][0] + tri[0][0]*tri[2][0] + tri[1][0]*tri[1][0] + tri[1][0]*tri[2][0] + tri[2][0]*tri[2][0])/6.;

  mat(2, 1) = area*(tri[0][1]*tri[0][1] + tri[0][1]*tri[1][1] + tri[0][1]*tri[2][1] + tri[1][1]*tri[1][1] + tri[1][1]*tri[2][1] + tri[2][1]*tri[2][1])/6.;

  return mat;

}
