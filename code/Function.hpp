#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <vector>
#include <initializer_list>
#include <stdexcept>
#include <cmath>

#include "ItgQuadForm.hpp"

namespace fct{
  typedef double R;
  typedef R (*RtoR)(R);
  typedef R (*R2toR)(R, R);


  void setepsilon(R);
  R getepsilon();

  class Function {
    const int p;//dimension of the domain
    const int q;//dimension of the codomain
    int deg = 2;//for integral calculus, the function will be seen as a polynomial of degree 'deg'

  public:
    Function() = delete;
    Function(int, int);
    Function(const Function&) = delete;
    Function(Function&&) = delete;
    virtual ~Function() = default;

    virtual std::vector<R> operator()(const std::vector<R>&) const = 0;
    std::vector<R> operator()(const Tri<2>&);
    Function& operator=(const Function&) = delete;
    Function& operator=(Function&&) = delete;
    virtual Function* clone() const = 0;
    const int getp() const;
    const int getq() const;
    void setdeg(int);

    R integrate(const Tri<2>&) const;
    ntuple<3, R> projection_on_P1(const Tri<2>&) const;
    R get_E_K(const Tri<2>&, ntuple<3, R>&) const;
    std::vector<Tri<2> > exportGnuplot(const std::vector<Tri<2> >&, const std::string& filename) const;
  };

  class LinCbnFct : public Function {
    Function** const fcts;
    R* const coeffs;
    const int nb;

  public:
    LinCbnFct() = delete;
    LinCbnFct(const LinCbnFct&);
    LinCbnFct(LinCbnFct&&) = delete;
    LinCbnFct(std::initializer_list<const Function*>, std::initializer_list<R>);
    ~LinCbnFct() override;

    std::vector<R> operator()(const std::vector<R>&) const override;
    LinCbnFct& operator=(const LinCbnFct&) = delete;
    LinCbnFct& operator=(LinCbnFct&&) = delete;
    LinCbnFct* clone() const override;

  };

  class CompFct : public Function {
    Function** const fcts;
    const int nb;

  public:
    CompFct() = delete;
    CompFct(const CompFct&);
    CompFct(CompFct&&) = delete;
    CompFct(std::initializer_list<const Function*>);
    ~CompFct() override;

    std::vector<R> operator()(const std::vector<R>&) const override;
    CompFct& operator=(const CompFct&) = delete;
    CompFct& operator=(CompFct&&) = delete;
    CompFct* clone() const override;

  };

  class RtoRCppFct : public Function {
    const RtoR fct;

  public:
    RtoRCppFct() = delete;
    RtoRCppFct(const RtoRCppFct&);
    RtoRCppFct(RtoRCppFct&&) = delete;
    RtoRCppFct(RtoR);
    ~RtoRCppFct() = default;

    std::vector<R> operator()(const std::vector<R>&) const override;
    RtoRCppFct& operator=(const RtoRCppFct&) = delete;
    RtoRCppFct& operator=(RtoRCppFct&&) = delete;
    RtoRCppFct* clone() const override;

  };

  class R2toRCppFct : public Function {
    const R2toR fct;

  public:
    R2toRCppFct() = delete;
    R2toRCppFct(const R2toRCppFct&);
    R2toRCppFct(R2toRCppFct&&) = delete;
    R2toRCppFct(R2toR);
    ~R2toRCppFct() = default;

    std::vector<R> operator()(const std::vector<R>&) const override;
    R2toRCppFct& operator=(const R2toRCppFct&) = delete;
    R2toRCppFct& operator=(R2toRCppFct&&) = delete;
    R2toRCppFct* clone() const override;

  };

  class ProdFct : public Function {
    Function** const fcts;
    const int nb;

  public:
    ProdFct() = delete;
    ProdFct(const ProdFct&);
    ProdFct(ProdFct&&) = delete;
    ProdFct(std::initializer_list<const Function*>);
    ~ProdFct();

    std::vector<R> operator()(const std::vector<R>&) const override;
    ProdFct& operator=(const ProdFct&) = delete;
    ProdFct& operator=(ProdFct&&) = delete;
    ProdFct* clone() const override;

  };

  class F_K : public Function {
    Tri<2> tri;

  public:
    F_K() = delete;
    F_K(const F_K&);
    F_K(F_K&&) = delete;
    F_K(const Tri<2>&);
    ~F_K() = default;

    std::vector<R> operator()(const std::vector<R>&) const override;
    F_K& operator=(const F_K&) = delete;
    F_K& operator=(F_K&&) = delete;
    F_K* clone() const override;

  };

  R square(R);
  R cube(R);
  R indR2(R, R);
  R p1(R, R);
  R p2(R, R);
  extern RtoRCppFct Square;
  extern RtoRCppFct Cube;
  extern R2toRCppFct IndR2;
  extern R2toRCppFct pr1;
  extern R2toRCppFct pr2;
}

#endif
