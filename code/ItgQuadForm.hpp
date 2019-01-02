#ifndef ITGQUADFORM_HPP
#define ITGQUADFORM_HPP

#include <vector>
#include <string>
#include <map>
#include <fstream>

#include "Tri.hpp"
#include "Matrix.hpp"


namespace fct{

  class ItgQuadForm{
  private:
    static std::map<int, std::vector<ntuple<3, double> > > pts_w;//points and weights for quadrature formulas
    static Matrix<3, 4> mat;//matrix for representing the linear system
    static const ItgQuadForm* const obj;

    ItgQuadForm();
    ItgQuadForm(const ItgQuadForm&) = delete;
    ItgQuadForm(ItgQuadForm&&) = delete;
    ~ItgQuadForm() = default;
    ItgQuadForm& operator=(const ItgQuadForm&) = delete;
    ItgQuadForm& operator=(ItgQuadForm&&) = delete;

  public:
    static const std::vector<ntuple<3, double> >& getpw(int);
    static Matrix<3, 4>& get_lin_sys(const Tri<2>&);
  };
}

#endif
