#include "Function.hpp"

using namespace std;
using namespace fct;

Function::Function(int pp, int qq) : p(pp), q(qq){}

const int Function::getp() const {
  return p;
}

const int Function::getq() const {
  return q;
}


LinCbnFct::LinCbnFct(const LinCbnFct& lcf) : Function(lcf.getp(), lcf.getq()), fcts(new Function*[lcf.nb]), coeffs(new R[lcf.nb]), nb(lcf.nb){

  for(int i=0;i<nb;i++){
    fcts[i] = lcf.fcts[i]->clone();
    coeffs[i] = lcf.coeffs[i];
  }

}

LinCbnFct::LinCbnFct(initializer_list<Function*> ilf, initializer_list<R> ilr) : Function(ilf.size()?ilf.begin()[0]->getp():(throw(logic_error("LinCbnFct({Function*}, {R}) : at least one function is expected\n"))), ilf.size()?ilf.begin()[0]->getq():0), fcts(new Function*[ilf.size()]), coeffs(new R[ilf.size()]), nb(ilf.size()){

  for(int i=1;i<nb;i++){
    if(getp()!=ilf.begin()[i]->getp() || getq()!=ilf.begin()[i]->getq()){
      delete fcts;
      delete coeffs;
      throw(logic_error("LinCbnFct({Function*}, {R}) : the functions must have the same domain and the same codomain\n"));
    }
  }

  for(int i=0;i<nb;i++){
    fcts[i] = ilf.begin()[i]->clone();
    if(i<ilr.size()){
      coeffs[i] = ilr.begin()[i];
    }else{
      coeffs[i] = (R)(1.0);
    }
  }

}

LinCbnFct::~LinCbnFct(){
  for(int i=0;i<nb;i++){
    delete fcts[i];
  }
  delete fcts;
  delete coeffs;
}

vector<R> LinCbnFct::operator()(const vector<R>& args) const{

  int values_nb = (getq()*args.size())/getp();
  vector<R> sums(values_nb, (R)(0.0));
  vector<R> aux;

  for(int i=0;i<nb;i++){
    aux = fcts[i]->operator()(args);
    for(int j=0;j<values_nb;j++){
      sums[j] += aux[j] * coeffs[i];
    }
  }

  return sums;
}

LinCbnFct* LinCbnFct::clone() const{
  return new LinCbnFct(*this);
}


CompFct::CompFct(const CompFct& cf) : Function(cf.getp(), cf.getq()), fcts(new Function*[cf.nb]), nb(cf.nb){
  for(int i=0;i<nb;i++){
    fcts[i] = cf.fcts[i]->clone();
  }
}

CompFct::CompFct(initializer_list<Function*> ilf) : Function(ilf.size()?ilf.begin()[ilf.size()-1]->getp():(throw(logic_error("CompFct({Function*}) : at least one function is expected\n")),1), ilf.size()?ilf.begin()[0]->getq():0), fcts(new Function*[ilf.size()]), nb(ilf.size()){

  for(int i=0;i<nb-1;i++){
    if(ilf.begin()[i]->getp() != ilf.begin()[i+1]->getq()){
      delete fcts;
      throw(logic_error("CompFct({Function*}) : composition is not consistent\n"));
    }
  }

  for(int i=0;i<nb;i++){
    fcts[nb-i-1] = ilf.begin()[i]->clone();
  }

}

CompFct::~CompFct(){

  for(int i=0;i<nb;i++){
    delete fcts[i];
  }
  delete fcts;

}

vector<R> CompFct::operator()(const vector<R>& args) const{
  vector<R> res(args.size());

  res = fcts[0]->operator()(args);
  for(int i=1;i<nb;i++){
    res = fcts[i]->operator()(res);
  }

  return res;

}

CompFct* CompFct::clone() const{
  return new CompFct(*this);
}


RtoRCppFct::RtoRCppFct(const RtoRCppFct& cppfct) : Function(cppfct.getp(), cppfct.getq()), fct(cppfct.fct){}

RtoRCppFct::RtoRCppFct(RtoR func) : Function(1, 1), fct(func){}

vector<R> RtoRCppFct::operator()(const vector<R>& args) const{

  vector<R> res(args.size());
  for(int i=0;i<args.size();i++){
    res[i] = (*fct)(args[i]);
  }

  return res;

}

RtoRCppFct* RtoRCppFct::clone() const{
  return new RtoRCppFct(*this);
}


R2toRCppFct::R2toRCppFct(const R2toRCppFct& cppfct) : Function(cppfct.getp(), cppfct.getq()), fct(cppfct.fct){}

R2toRCppFct::R2toRCppFct(R2toR func) : Function(2, 1), fct(func){}

vector<R> R2toRCppFct::operator()(const vector<R>& args) const{

  if(args.size()%2 != 0){
    throw(invalid_argument("R2toRCppFct::operator() : the argument's size must be even\n"));
  }

  int out_nb = args.size() >> 1;
  vector<R> res(out_nb);
  for(int i=0;i<out_nb;i++){
    res[i] = (*fct)(args[i<<1], args[(i<<1)+1]);
  }

  return res;

}

R2toRCppFct* R2toRCppFct::clone() const{
  return new R2toRCppFct(*this);
}
