#ifndef NTUPLE_HPP
#define NTUPLE_HPP

#include <stdexcept>
#include <iostream>
#include <initializer_list>

template<int n, typename NUM>
class ntuple;

template<int n, typename NUM>
ntuple<n, NUM> operator+(const ntuple<n, NUM>&, const ntuple<n, NUM>&);

template<int n, typename NUM>
ntuple<n, NUM> operator*(const ntuple<n, NUM>&, const ntuple<n, NUM>&);

template<int n, typename NUM>
ntuple<n, NUM> operator*(const NUM&, const ntuple<n, NUM>&);

template<int n, typename NUM>
std::ostream& operator<<(std::ostream&, const ntuple<n, NUM>&);



template<int n, typename NUM>
class ntuple {

  NUM t[n];

public:

  ntuple()=default;
  ntuple(const ntuple&);
  ntuple(ntuple&&);
  ntuple(std::initializer_list<NUM>);
  ~ntuple()=default;
  
  ntuple& operator=(const ntuple&);
  ntuple& operator=(ntuple&&);
  NUM& operator[] (int);
  const NUM& operator[] (int) const;
  friend ntuple operator+<>(const ntuple&, const ntuple&);
  friend ntuple operator*<>(const ntuple&, const ntuple&);
  friend ntuple operator*<>(const NUM&, const ntuple&);
  friend std::ostream& operator<<<>(std::ostream&, const ntuple&);
};

template<int n, typename NUM>
ntuple<n, NUM>::ntuple(const ntuple<n, NUM>& nt){
  for(int i=0;i<n;i++){
    t[i] = nt.t[i];
  }
}

template<int n, typename NUM>
ntuple<n, NUM>::ntuple(ntuple<n, NUM>&& nt){
  for(int i=0;i<n;i++){
    t[i] = nt.t[i];
  }
}

template<int n, typename NUM>
ntuple<n, NUM>::ntuple(std::initializer_list<NUM> il){
  if(il.size()>n){
    throw(std::invalid_argument("ntuple({NUM}) : too many elements in the list\n"));
  }
  for(int i=0;i<il.size();i++){
    t[i] = il.begin()[i];
  }
}

template<int n, typename NUM>
ntuple<n, NUM>& ntuple<n, NUM>::operator=(const ntuple<n, NUM>& nt){
  for(int i=0;i<n;i++){
    t[i] = nt.t[i];
  }
  return *this;
}

template<int n, typename NUM>
ntuple<n, NUM>& ntuple<n, NUM>::operator=(ntuple<n, NUM>&& nt){
  for(int i=0;i<n;i++){
    t[i] = nt.t[i];
  }
  return *this;
}

template<int n, typename NUM>
NUM& ntuple<n, NUM>::operator[](int i){
  if(i<0 || i>=n){
    throw std::out_of_range("bad subscript for class ntuple");
  }
  return t[i];
}

template<int n, typename NUM>
const NUM& ntuple<n, NUM>::operator[](int i) const{
  if(i<0 || i>=n){
    throw std::out_of_range("bad subscript for class ntuple");
  }
  return t[i];
}

template<int n, typename NUM>
ntuple<n, NUM> operator+(const ntuple<n, NUM>& op1, const ntuple<n, NUM>& op2){
  ntuple<n, NUM> res;
  for(int i=0;i<n;i++){
    res.t[i] = op1.t[i] + op2.t[i];
  }
  return res;
}

template<int n, typename NUM>
ntuple<n, NUM> operator*(const ntuple<n, NUM>& op1, const ntuple<n, NUM>& op2){
  ntuple<n, NUM> res;
  for(int i=0;i<n;i++){
    res.t[i] = op1.t[i] * op2.t[i];
  }
  return res;
}

template<int n, typename NUM>
ntuple<n, NUM> operator*(const NUM& op1, const ntuple<n, NUM>& op2){
  ntuple<n, NUM> res;
  for(int i=0;i<n;i++){
    res.t[i] = op1 * op2.t[i];
  }
  return res;
}

template<int n, typename NUM>
std::ostream& operator<<(std::ostream& s, const ntuple<n, NUM>& nt){
  s << nt.t[0];
  for(int i=1;i<n;i++){
    s << " " << nt.t[i];
  }
  return s;
}

#endif
