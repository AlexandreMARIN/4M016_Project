#include<iostream>

#include "../code/Matrix.hpp"

using namespace std;

int main(int argc, char* argv[]){

  Matrix<3, 4> mat;
  vector<ntuple<2, int> > ldgent;

  mat(0, 0) = 2.;
  mat(0, 1) = 4.;
  mat(1, 1) = 2.;
  mat(1, 0) = 2.5;
  mat(2, 2) = 0.5;
  mat(0, 3) = 1;
  mat(1, 3) = 2.;
  mat(2, 3) = 6.;
  mat(2, 1) = 2.;

  cout << "mat =\n" << mat << "\n";
  cout << "row-reduced echelon form :\n" << mat.gauss(ldgent);

  cout << "\nleading entries for gaussian elimination :\n";
  for(int i=0;i<ldgent.size();i++){
    cout << "(" << ldgent[i] << ")\n";
  }

  cout << "\n\n";

  return 0;
}
