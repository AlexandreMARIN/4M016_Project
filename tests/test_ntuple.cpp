#include <iostream>

#include "../code/ntuple.hpp"

using namespace std;

template ostream& operator<<(ostream&, const ntuple<3, double>&);

int main(int argc, char* argv[]){

  ntuple<3, double>a{0.5, -0.5, 1.5}, b{0., 1., 2.5};
  const ntuple<3, double>c(a);

  cout << "a = (" << a << ")\nb = (" << b << ")\n";

  cout << "c = (" << c << ")\nc[1] = " << c[1] << "\n";

  cout << "a+b = (" << a+b << ")\na*b = (" << a*b << ")\n";

  cout << "10*a = (" << 10.*a << ")\n\n";

  return 0;

}
