#include "Function.hpp"

#include <queue>
#include <list>
#include <fstream>

using namespace std;
using namespace fct;

namespace fct{
  static R epsilon = 1e-9;
}
constexpr R min_eps = 1e-10;


void fct::setepsilon(R e){
  if(e>min_eps && e < 1.0){
    epsilon = e;
  }
}

R fct::getepsilon(){
  return epsilon;
}

Function::Function(int pp, int qq) : p(pp), q(qq){}

vector<R> Function::operator()(const Tri<2>& tri){
  vector<R> args{tri[0][0], tri[0][1], tri[1][0], tri[1][1], tri[2][0], tri[2][1]};
  return operator()(args);
}

const int Function::getp() const {
  return p;
}

const int Function::getq() const {
  return q;
}

void Function::setdeg(int newdeg){
  if(newdeg>0){
    deg = newdeg;
  }
}

R Function::integrate(const Tri<2>& tri) const{
  if(getp()!=2 || getq()!=1){
    throw(logic_error("Function::integrate(const Tri<2>&) : the implicit argument must be a function of that sort : R^2 -> R\n"));
  }

  const vector<ntuple<3, double> >& pw = ItgQuadForm::getpw(deg);
  F_K f_k(tri);

  CompFct aux({this, &f_k});
  int npts = pw.size();
  vector<R> args(npts<<1);
  vector<R> res;
  R integ;
  R area = 0.5*abs((tri[1][0]-tri[0][0])*(tri[2][1]-tri[0][1]) - (tri[1][1]-tri[0][1])*(tri[2][0]-tri[0][0]));

  for(int i=0;i<npts;i++){
    args[i<<1] = pw[i][0];
    args[(i<<1)+1] = pw[i][1];
  }

  res = aux(args);

  integ = res[0] * pw[0][2];
  for(int i=1;i<npts;i++){
    integ += res[i] * pw[i][2];
  }

  return area*((R)0.5)*integ;
}

ntuple<3, R> Function::projection_on_P1(const Tri<2>& tri) const{
  Matrix<3, 4>& sys = ItgQuadForm::get_lin_sys(tri);
  ntuple<3, R> proj;
  vector<ntuple<2, int> > ldgent;
  sys(0, 3) = integrate(tri);

  ProdFct fx({this, &pr1}), fy({this, &pr2});
  fx.setdeg(deg+1);
  fy.setdeg(deg+1);
  sys(1, 3) = fx.integrate(tri);
  sys(2, 3) = fy.integrate(tri);

  Matrix<3, 4> copy = sys;//for the exception below

  sys.gauss(ldgent);
  if(ldgent.size()<3 || ldgent[0][1]==3 || ldgent[1][1]==3 || ldgent[2][1]==3){
    cerr << "\nerror in ntuple<3, R> Function::projection_on_P1(const Tri<2>& tri):\n\nprojection on triangle:\n{( " << tri[0] << " ), ( " << tri[1] << " ), ( " << tri[2] << " )}\n\nmatrix before and after Gauss algorithm:\n" << copy << "~\n" << sys << "\n\n";
    throw(range_error("Function::projection_on_P1 : the matrix is not invertible\n"));
  }

  proj[ldgent[0][1]] = sys(ldgent[0][0], 3);
  proj[ldgent[1][1]] = sys(ldgent[1][0], 3);
  proj[ldgent[2][1]] = sys(ldgent[2][0], 3);

  return proj;
}

R Function::get_E_K(const Tri<2>& tri, ntuple<3, R>& proj) const{
  proj = projection_on_P1(tri);
  LinCbnFct diff({this, &pr1, &pr2, &IndR2}, {-1.0, proj[0], proj[1], proj[2]});
  CompFct err({&Square, &diff});

  return err.integrate(tri);
}

vector<Tri<2> > Function::exportGnuplot(const vector<Tri<2> >& mesh, const string& filename) const{
  ofstream file(filename);
  vector<Tri<2> > newmesh((int)(mesh.size()*1.5));
  queue<Tri<2>, list<Tri<2> > > qt;
  Tri<2> curtri;
  Tri<2> splittri[6];
  int option, maxlg = 10000, lg = 0;
  R errors[6];
  R total_err[3];
  vector<R> res(3);
  ntuple<3, R> proj[6];

  //insert triangles in the queue if the resulting errors are greater than epsilon
  for(const Tri<2>& tri : mesh){
    if(get_E_K(tri, proj[0])<epsilon){
      res[0] = proj[0][0]*tri[0][0] + proj[0][1]*tri[0][1] + proj[0][2];
      res[1] = proj[0][0]*tri[1][0] + proj[0][1]*tri[1][1] + proj[0][2];
      res[2] = proj[0][0]*tri[2][0] + proj[0][1]*tri[2][1] + proj[0][2];
      file << Tri<3>(tri, res) << "\n\n";
      lg++;
      newmesh.push_back(tri);
    }else{
      qt.push(tri);
    }
  }

  //here we refine the mesh
  while(!qt.empty()){
    if(lg>maxlg){
      throw(runtime_error("Function::exportGnuplot() : the file "+filename+" is too big\n"));
    }
    curtri = qt.front();
    qt.pop();
    //we split the current triangle
    curtri.split(splittri);
    for(int i=0;i<6;i++){
      errors[i] = get_E_K(splittri[i], proj[i]);
    }
    total_err[0] = errors[0] + errors[1];
    total_err[1] = errors[2] + errors[3];
    total_err[2] = errors[4] + errors[5];
    //we choose the split which minimizes as much as possible the total error
    if(total_err[0]<total_err[1]){
      if(total_err[0]<total_err[2]){
	option = 0;
      }else{
	option = 2;
      }
    }else{
      if(total_err[1]<total_err[2]){
	option = 1;
      }else{
	option = 2;
      }
    }

    //then we insert into the queue these triangles if error E_K is greater than epsilon
    option <<= 1;
    if(errors[option]<epsilon){
      res[0] = proj[option][0]*splittri[option][0][0] + proj[option][1]*splittri[option][0][1] + proj[option][2];
      res[1] = proj[option][0]*splittri[option][1][0] + proj[option][1]*splittri[option][1][1] + proj[option][2];
      res[2] = proj[option][0]*splittri[option][2][0] + proj[option][1]*splittri[option][2][1] + proj[option][2];
      file << Tri<3>(splittri[option], res) << "\n\n";
      lg++;
      newmesh.push_back(splittri[option]);
    }else{
      qt.push(splittri[option]);
    }
    option++;

    if(errors[option]<epsilon){
      res[0] = proj[option][0]*splittri[option][0][0] + proj[option][1]*splittri[option][0][1] + proj[option][2];
      res[1] = proj[option][0]*splittri[option][1][0] + proj[option][1]*splittri[option][1][1] + proj[option][2];
      res[2] = proj[option][0]*splittri[option][2][0] + proj[option][1]*splittri[option][2][1] + proj[option][2];
      file << Tri<3>(splittri[option], res) << "\n\n";
      lg++;
      newmesh.push_back(splittri[option]);
    }else{
      qt.push(splittri[option]);
    }
  }

  file.close();
  return newmesh;
}


LinCbnFct::LinCbnFct(const LinCbnFct& lcf) : Function(lcf.getp(), lcf.getq()), fcts(new Function*[lcf.nb]), coeffs(new R[lcf.nb]), nb(lcf.nb){

  for(int i=0;i<nb;i++){
    fcts[i] = lcf.fcts[i]->clone();
    coeffs[i] = lcf.coeffs[i];
  }

}

LinCbnFct::LinCbnFct(initializer_list<const Function*> ilf, initializer_list<R> ilr) : Function(ilf.size()?ilf.begin()[0]->getp():(throw(logic_error("LinCbnFct({Function*}, {R}) : at least one function is expected\n"))), ilf.size()?ilf.begin()[0]->getq():0), fcts(new Function*[ilf.size()]), coeffs(new R[ilf.size()]), nb(ilf.size()){

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

CompFct::CompFct(initializer_list<const Function*> ilf) : Function(ilf.size()?ilf.begin()[ilf.size()-1]->getp():(throw(logic_error("CompFct({Function*}) : at least one function is expected\n")),1), ilf.size()?ilf.begin()[0]->getq():0), fcts(new Function*[ilf.size()]), nb(ilf.size()){

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


ProdFct::ProdFct(const ProdFct& pf) : Function(pf.getp(), pf.getq()), fcts(new Function*[pf.nb]), nb(pf.nb){
  for(int i=0;i<nb;i++){
    fcts[i] = pf.fcts[i]->clone();
  }
}

ProdFct::ProdFct(initializer_list<const Function*> ilf) : Function(ilf.size()!=0?ilf.begin()[0]->getp():-1, 1), fcts(new Function*[ilf.size()]), nb(ilf.size()){

  if(!ilf.size()){
    throw(invalid_argument("ProdFct({Function*}) : at least one function is expected\n"));
  }
  for(int i=1;i<nb;i++){
    if(ilf.begin()[i]->getp()!=getp()){
      delete fcts;
      throw(logic_error("ProdFct({Function*}) : the functions must have the same domain\n"));
    }
  }
  for(int i=0;i<nb;i++){
    if(ilf.begin()[i]->getq()!=1){
      delete fcts;
      throw(logic_error("ProdFct({Function*}) : the functions must have R for codomain (i.e. getq()==1)\n"));
    }
  }

  for(int i=0;i<nb;i++){
    fcts[i] = ilf.begin()[i]->clone();
  }

}

ProdFct::~ProdFct(){
  for(int i=0;i<nb;i++){
    delete fcts[i];
  }
  delete fcts;
}

vector<R> ProdFct::operator()(const vector<R>& args) const{

  vector<R> res(fcts[0]->operator()(args));
  vector<R> aux;

  for(int i=1;i<nb;i++){
    aux = fcts[i]->operator()(args);
    for(int j=0;j<res.size();j++){
      res[j] *= aux[j];
    }
  }

  return res;
}

ProdFct* ProdFct::clone() const{
  return new ProdFct(*this);
}


F_K::F_K(const F_K& f) : Function(2, 2), tri(f.tri){}

F_K::F_K(const Tri<2>& t) : Function(2, 2), tri(t){}

vector<R> F_K::operator()(const vector<R>& args) const{

  if(args.size()%2 != 0){
    throw(invalid_argument("F_K::operator() : argument's size must be even\n"));
  }

  vector<R> res(args.size());

  for(int i=0;i<res.size();i+=2){
    res[i] = (((R)1.0)-args[i]-args[i+1])*tri[0][0] + args[i]*tri[1][0] + args[i+1]*tri[2][0];
    res[i+1] = (((R)1.0)-args[i]-args[i+1])*tri[0][1] + args[i]*tri[1][1] + args[i+1]*tri[2][1];
  }

  return res;
}

F_K* F_K::clone() const{
  return new F_K(*this);
}



R fct::square(R x){
  return x*x;
}

R fct::cube(R x){
  return x*x*x;
}

R fct::indR2(R x, R y){
  return 1.0;
}

R fct::p1(R x, R y){
  return x;
}

R fct::p2(R x, R y){
  return y;
}

RtoRCppFct fct::Square(&square);
RtoRCppFct fct::Cube(&cube);
R2toRCppFct fct::IndR2(&indR2);
R2toRCppFct fct::pr1(&p1);
R2toRCppFct fct::pr2(&p2);
