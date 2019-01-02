#ifndef TRI_HPP
#define TRI_HPP

#include "ntuple.hpp"

template<int dim> class Tri;
template<int dim>
std::ostream& operator<<(std::ostream&, Tri<dim>&);

template<int dim>
class Tri{
  static_assert(dim==2 ||dim==3, "Tri : dim must be 2 or 3\n");
  static int NbTri;
  ntuple<dim, double> vert[3];

public:
  //constructors
  Tri();
  Tri(const Tri&);
  Tri(Tri&&);
  Tri(ntuple<dim, double>, ntuple<dim, double>, ntuple<dim, double>);

  ~Tri();

  //operators
  Tri& operator=(const Tri&);
  Tri& operator=(Tri&&);
  friend std::ostream& operator<<<>(std::ostream&, Tri&);
  const ntuple<dim, double>& operator[](int) const;
  ntuple<dim, double>& operator[](int);

  //methods
  static int Nb();
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
std::ostream& operator<<(std::ostream& os, Tri<dim>& tri){
  os << "{(" << tri.vert[0] << "),\n (" << tri.vert[1] << "),\n (" << tri.vert[2] << ")}\n";
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
