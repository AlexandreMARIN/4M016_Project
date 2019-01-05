#include "Function.hpp"

using namespace std;
using namespace fct;

int main(int argc, char* argv[]){

  string arg1, arg2;
  int nbt = 4;
  vector<Tri<2> > mesh;

  try{
    if(argc>1){
      arg1 = argv[1];
      setepsilon(stod(arg1, nullptr));
      if(argc>2){
	arg2 = argv[2];
	nbt = stoi(arg2, nullptr);
      }
    }
  }

  catch(out_of_range&){
    cerr << "You have given a value which is out of range for \"double\"\n";
    return 1;
  }
  catch(invalid_argument&){
    cerr << "the given string cannot be converted to a \"double\"\n";
    return 1;
  }

  cout << "\nThis program will run with:\n- a mesh of [-1, 1]^2 composed of " << nbt << " triangles\n- epsilon = " << getepsilon() << "\n";

  //here we build the functions f_1 and f_2
  CompFct x2({&Square, &pr1}), y2({&Square, &pr2}), y3({&Cube, &pr2});
  LinCbnFct f_1({&x2, &y2}, {1.0, 1.0});

  LinCbnFct func({&pr1, &pr2}, {2.0, 2.0});
  RtoRCppFct Sin(&sin), Tanh(&tanh);
  LinCbnFct func2({&Sin}, {5.0});
  CompFct func3({&Tanh, &func2, &func});
  LinCbnFct f_2({&x2, &y3, &func3}, {});

  mesh = Tri<2>::generateMesh(-1., 1., -1., 1., nbt);
  f_2.setdeg(5);

  f_1.exportGnuplot(mesh, "f_1.data");
  f_2.exportGnuplot(mesh, "f_2.data");

  return 0;
}
