#include "Tri.hpp"

using namespace std;

template class Tri<2>;
template class Tri<3>;

template<int dim> int Tri<dim>::NbTri = 0;

template<>
template<>
Tri<3>::Tri(const Tri<2>& tri) : Tri(){
  vert[0] = ntuple<3, double>({tri[0][0], tri[0][1], 0.0});
  vert[1] = ntuple<3, double>({tri[1][0], tri[1][1], 0.0});
  vert[2] = ntuple<3, double>({tri[2][0], tri[2][1], 0.0});
}

template<>
template<>
Tri<3>::Tri(const Tri<2>& tri2, const std::vector<double>& v) : Tri(){
  if(v.size()!=3){
    throw(logic_error("Tri<3>::Tri(const Tri<2>& tri2, const std::vector<double>& v) : v's size must be 3\n"));
  }

  vert[0] = ntuple<3, double>({tri2[0][0], tri2[0][1], v[0]});
  vert[1] = ntuple<3, double>({tri2[1][0], tri2[1][1], v[1]});
  vert[2] = ntuple<3, double>({tri2[2][0], tri2[2][1], v[2]});
}

template<>
void Tri<2>::split(Tri<2> tris[6]) const{
  for(int i=0;i<3;i++){
    tris[(i<<1)].vert[i] = vert[i];
    tris[(i<<1)].vert[(i+1)%3] = vert[(i+1)%3];
    tris[(i<<1)].vert[(i+2)%3] = 0.5*(vert[(i+1)%3]+vert[(i+2)%3]);

    tris[(i<<1)+1].vert[i] = vert[i];
    tris[(i<<1)+1].vert[(i+1)%3] = vert[(i+2)%3];
    tris[(i<<1)+1].vert[(i+2)%3] = tris[(i<<1)].vert[(i+2)%3];
  }
}

template<>
vector<Tri<2> > Tri<2>::generateMesh(double a, double b, double c, double d, int nb){

  if(a>b || c>d || nb<=0){
    throw(invalid_argument("vector<Tri<2> > Tri<2>::generateMesh(double a, double b, double c, double d, int nb) : we must have: a<=b, c<=d and nb>0\n"));
  }

  ntuple<2, double> pM({a, d}), pN({b, d}), pO({b, c}), pP({a, c}), pQ({0.5*(a+b), 0.5*(c+d)});
  vector<Tri<2> > mesh{Tri(pQ, pM, pN), Tri(pQ, pN, pO), Tri(pQ, pO, pP), Tri(pQ, pP, pM)};
  vector<Tri<2> > aux(nb);
  Tri<2> ntri;

  mesh.reserve((nb<4)?4:nb);

  while(mesh.size()<nb){
    aux.resize(0);
    for(Tri<2>& tri : mesh){
      ntri.vert[0] = 0.5*(tri.vert[1]+tri.vert[2]);
      ntri.vert[1] = tri.vert[1];
      ntri.vert[2] = tri.vert[0];
      aux.push_back(ntri);
      ntri.vert[1] = tri.vert[2];
      aux.push_back(ntri);
    }
    mesh.swap(aux);
  }

  return mesh;
}
