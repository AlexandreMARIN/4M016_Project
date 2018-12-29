#include <iostream>
#include <cstdlib>
#include <cmath>

#include "../code/Function.hpp"

using namespace std;
using namespace fct;

R mean(R x, R y){
  return 0.5*(x+y);
}

int main(int argc, char* argv[]){

  int option = 1;

  R2toRCppFct Mean(&mean);
  RtoRCppFct Exp(&exp);
  vector<R> xargs{0., 2., 1.};

  if(argc>1){
    option = atoi(argv[1]);
  }

  try{
    switch(option){
    case 1:
      LinCbnFct({&Exp, &Mean}, {});
      break;
    case 2:
      CompFct({&Mean, &Exp});
      break;
    case 3:
      LinCbnFct({}, {1});
      break;
    case 4:
      CompFct({});
      break;
    default:
      Mean(xargs);
    }
  }

  catch(invalid_argument& ia){
    cerr << ia.what();
  }

  catch(logic_error& le){
    cerr << le.what();
  }

  return 0;
}
