#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <set>

#include "ntuple.hpp"

template<int m, int n>
class Matrix;

template<int m, int n>
std::ostream& operator<<(std::ostream&, const Matrix<m, n>&);

template<int m, int n>
class Matrix{
  static_assert(m<=10 && n<=10, "Matrix : integers are too large\n");
  static_assert(m>0 && n>0, "Matrix : integers must be positive\n");
  static_assert(m!=1 || n!=1, "Matrix : 1 by 1 matrices are forbidden\n");

  double ent[m][n];//entries

public:
  Matrix();
  Matrix(const Matrix&);
  Matrix(Matrix&&) = delete;
  ~Matrix() = default;

  Matrix& operator=(const Matrix&);
  Matrix& operator=(Matrix&&) = delete;

  double& operator()(int, int);
  const double& operator()(int, int) const;
  friend std::ostream& operator<<<>(std::ostream&, const Matrix&);

  Matrix& gauss(std::vector<ntuple<2, int> >&);
};

template<int m, int n>
Matrix<m, n>::Matrix(){
  int i, j;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      ent[i][j] = 0.0;
    }
  }
}

template<int m, int n>
Matrix<m, n>::Matrix(const Matrix<m, n>& mat){
  int i, j;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      ent[i][j] = mat.ent[i][j];
    }
  }
}

template<int m, int n>
Matrix<m, n>& Matrix<m, n>::operator=(const Matrix<m, n>& mat){
  int i, j;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      ent[i][j] = mat.ent[i][j];
    }
  }
}

template<int m, int n>
double& Matrix<m, n>::operator()(int i, int j){
  if(i<0 || i>=m){
    throw(std::out_of_range("Matrix::operator(int, int) : row index must be in [0, m-1]\n"));
  }
  if(j<0 || j>=n){
    throw(std::out_of_range("Matrix::operator(int, int) : column index must be in [0, n-1]\n"));
  }

  return ent[i][j];
}

template<int m, int n>
const double& Matrix<m, n>::operator()(int i, int j) const{
  if(i<0 || i>=m){
    throw(std::out_of_range("Matrix::operator(int, int) : row index must be in [0, m-1]\n"));
  }
  if(j<0 || j>=n){
    throw(std::out_of_range("Matrix::operator(int, int) : column index must be in [0, n-1]\n"));
  }

  return ent[i][j];
}

template<int m, int n>
std::ostream& operator<<(std::ostream& os, const Matrix<m, n>& mat){
  int i, j;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      os << mat(i, j) << " ";
    }
    os << "\n";
  }
  os << "\n";
}

template<int m, int n>
Matrix<m, n>& Matrix<m, n>::gauss(std::vector<ntuple<2, int> >& ldgent){

  //ldgent will contain positions of the leading entries
  ldgent.resize((m<n)?m:n);
  std::set<int> row_ind;
  int i, j, k, l, ldg1, ldg2;
  int ldgnb = 0;

  //at ther beginning of this method, row_ind must contain indices of non-null rows
  //at the end of this method, row_ind will contain indices of rows which do not contain any leading entry
  for(i=0;i<m;i++){
    row_ind.insert(i);
  }

  //we suppress indices of null rows
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      if(ent[i][j]!=0.){
	break;
      }
    }
    if(j==n){
      row_ind.erase(i);
    }
  }

  //for each column
  for(j=0;j<n;j++){
    i = -1;
    //for each row which has not been used for pivoting
    for(auto ind : row_ind){
      if(ent[ind][j]!=0.){
	i = ind;
	break;
      }
    }
    if(i>=0){
      //we are going to use the row i for pivoting
      row_ind.erase(i);
      for(l=j+1;l<n;l++){
	ent[i][l] /= ent[i][j];
      }
      ent[i][j] = 1.;
      for(auto ind2 : row_ind){
	for(l=j+1;l<n;l++){
	  ent[ind2][l] -= ent[ind2][j]*ent[i][l];
	}
	ent[ind2][j] = 0.;
      }
      ldgent[ldgnb++] = ntuple<2, int>({i, j});
    }
  }

  //some tranvections to put zeros above leading entries
  for(ldg1=1;ldg1<ldgnb;ldg1++){
    for(ldg2=0;ldg2<ldg1;ldg2++){
      i = ldgent[ldg1][0];
      j = ldgent[ldg1][1];
      k = ldgent[ldg2][0];
      for(l=j+1;l<n;l++){
	ent[k][l] -= ent[k][j]*ent[i][l];
      }
      ent[k][j] = 0.;
    }
  }

  ldgent.resize(ldgnb);
  return *this;
}

#endif
