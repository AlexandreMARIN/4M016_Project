#include "Tri.hpp"

using namespace std;

template class Tri<2>;
template class Tri<3>;

template<int dim> int Tri<dim>::NbTri = 0;
