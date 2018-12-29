#include <iostream>
#include <cstdlib>
#include <cmath>

#include "../code/Function.hpp"

using namespace std;
using namespace fct;

R sq(R x){
  return x*x;
}

int main(int argc, char* argv[]){

  vector<R> xargs(argc-1);

  RtoRCppFct Cos{&cos};
  RtoRCppFct Sin{&sin};
  RtoRCppFct Sq{&sq};

  CompFct Cos2({&Sq, &Cos});
  CompFct Sin2({&Sq, &Sin});
  LinCbnFct One({&Cos2, &Sin2}, {});

  if(argc>1){
    for(int i=1;i<argc;i++){
      xargs[i-1] = atof(argv[i]);
    }
  }else{
    xargs = vector<R>{0, 1, M_PI};
  }


  vector<R> res = Cos(xargs);


  for(int i=0;i<xargs.size();i++){
    cout << "cos(" << xargs[i] << ") = " << res[i] << "\n";
  }

  res = One(xargs);

  cout << "\nwith f(x) = cos^2(x)+sin^2(x) = 1,\n";
  for(int i=0;i<xargs.size();i++){
    cout << "f(" << xargs[i] << ") = " << res[i] << "\n";
  }



  return 0;
}
