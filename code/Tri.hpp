#ifndef TRI_HPP
#define TRI_HPP

#include <vector>
#include "ntuple.hpp"

template<int dim> class Tri;
template<int dim>
std::ostream& operator<<(std::ostream&, const Tri<dim>&);


template<int dim>
class Tri{
  static_assert(dim==2 ||dim==3, "Tri : dim must be 2 or 3\n");
  static int NbTri;
  ntuple<dim, double> vert[3];//vertices

public:
  //constructors
  Tri();
  Tri(const Tri&);
  Tri(Tri&&);
  Tri(ntuple<dim, double>, ntuple<dim, double>, ntuple<dim, double>);
  template<int dim2> Tri(const Tri<dim2>&);//it will be specialized
  template<int dim2> Tri(const Tri<dim2>&, const std::vector<double>&);//it will be specialized

  ~Tri();

  //operators
  Tri& operator=(const Tri&);
  Tri& operator=(Tri&&);
  friend std::ostream& operator<<<>(std::ostream&, const Tri&);
  const ntuple<dim, double>& operator[](int) const;
  ntuple<dim, double>& operator[](int);

  //methods
  static int Nb();
  void split(Tri[6]) const;//it will be specialized

  static std::vector<Tri> generateMesh(double, double, double, double, int);//it will be specialized

};



template<int dim>
Tri<dim>::Tri(){
  NbTri++;
}

template<int dim>
Tri<dim>::Tri(const Tri<dim>& tri){
  vert[0] = tri.vert[0];
  vert[1] = tri.vert[1];
  vert[2] = tri.vert[2];
  NbTri++;
}

template<int dim>
Tri<dim>::Tri(Tri<dim>&& tri){
  vert[0] = tri.vert[0];
  vert[1] = tri.vert[1];
  vert[2] = tri.vert[2];
  NbTri++;
}

template<int dim>
Tri<dim>::Tri(ntuple<dim, double> a, ntuple<dim, double> b, ntuple<dim, double> c){
  vert[0] = a;
  vert[1] = b;
  vert[2] = c;
  NbTri++;
}

template<int dim>
Tri<dim>::~Tri(){
  NbTri--;
}

template<int dim>
Tri<dim>& Tri<dim>::operator=(const Tri<dim>& tri){
  vert[0] = tri.vert[0];
  vert[1] = tri.vert[1];
  vert[2] = tri.vert[2];
  return *this;
}

template<int dim>
Tri<dim>& Tri<dim>::operator=(Tri<dim>&& tri){
  vert[0] = tri.vert[0];
  vert[1] = tri.vert[1];
  vert[2] = tri.vert[2];
  return *this;
}

template<int dim>
std::ostream& operator<<(std::ostream& os, const Tri<dim>& tri){
  //os << "{(" << tri.vert[0] << "),\n (" << tri.vert[1] << "),\n (" << tri.vert[2] << ")}\n";
  //above, there is a test version
  os << tri.vert[0] << "\n" << tri.vert[1] << "\n" << tri.vert[2] << "\n" << tri.vert[0] << "\n";
  return os;
}

template<int dim>
inline const ntuple<dim, double>& Tri<dim>::operator[](int i) const{

  switch(i){
  case 0:
  case 1:
  case 2:
    return vert[i];
    break;
  default:
    throw std::out_of_range("index not in [0, 2] in operator[] for class Tri");
  }
}

template<int dim>
inline ntuple<dim, double>& Tri<dim>::operator[](int i){

  switch(i){
  case 0:
  case 1:
  case 2:
    return vert[i];
    break;
  default:
    throw std::out_of_range("index not in [0, 2] in operator[] for class Tri");
  }
}

template<int dim>
int Tri<dim>::Nb(){
  return NbTri;
}



#endif
