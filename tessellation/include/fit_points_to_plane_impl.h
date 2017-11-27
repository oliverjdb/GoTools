#ifndef _FIT_POINTS_TO_PLANE_IMPL_H
#define _FIT_POINTS_TO_PLANE_IMPL_H

#include <cmath>
#include <algorithm>
#include <vector>
#include "GoTools/utils/GeneralFunctionMinimizer.h"

namespace TesselateUtils {

// ----------------------------------------------------------------------------  
template<typename P>
std::array<double, 4> fit_points_to_plane(const P* const points,
                                          const unsigned int num_points)
// ----------------------------------------------------------------------------  
{
  // @@ Since we do not have a linear algebra package easily available, we do this the
  // roundabout way - we find the angles (theta, psi) such that the plane having the normal
  // [cos(theta) cos(psi), sin(theta) cos(psi), sin(psi)] minimize the squared distance to the
  // points.  We use a function minimizer for this purpose.

  // First, we translate the points to the origin
  P mean_pt = std::accumulate(points, points+num_points, P {0.0, 0.0, 0.0});
  mean_pt /= (double) num_points;

  std::vector<P> cpts(points, points + num_points); // centered points
  std::transform(cpts.begin(), cpts.end(), cpts.begin(),
                 [mean_pt](const P& p) {return p - mean_pt;});

  
  const double a11 = std::accumulate(cpts.begin(), cpts.end(), 0.0,
                                     [](double acc, const P& p) {return acc + p[0] * p[0];});
  const double a12 = std::accumulate(cpts.begin(), cpts.end(), 0.0,
                                     [](double acc, const P& p) {return acc + p[0] * p[1];});
  const double a13 = std::accumulate(cpts.begin(), cpts.end(), 0.0,
                                     [](double acc, const P& p) {return acc + p[0] * p[2];});
  const double a22 = std::accumulate(cpts.begin(), cpts.end(), 0.0,
                                     [](double acc, const P& p) {return acc + p[1] * p[1];});
  const double a23 = std::accumulate(cpts.begin(), cpts.end(), 0.0,
                                     [](double acc, const P& p) {return acc + p[1] * p[2];});
  const double a33 = std::accumulate(cpts.begin(), cpts.end(), 0.0,
                                     [](double acc, const P& p) {return acc + p[2] * p[2];});

  P A1 {a11, a12, a13}; // first row (or column) of symmetric matrix A
  P A2 {a12, a22, a23}; // second row (or column)
  P A3 {a13, a23, a33}; // third row (or column)

  // our function to minimize is f(theta, psi) = p(theta,psi)' A P(theta,psi)
  
  struct TmpFunc {
    
    P A1, A2, A3;
    TmpFunc(P a1, P a2, P a3) : A1(a1), A2(a2), A3(a3) {}
    double operator()(const double* arg) const {
      // compute function value
      const P n = nvec(arg[0], arg[1]);
      return n[0] * scalprod(A1, n) + n[1] * scalprod(A2, n) +  n[2] * scalprod(A3, n);
    }
    void grad(const double* arg, double* grad) const {
      const P n = nvec(arg[0], arg[1]);      
      const P nt = nvec_dtheta(arg[0], arg[1]);
      const P np = nvec_dpsi(arg[0], arg[1]);
      
      grad[0] = 2 * (nt[0] * scalprod(A1, n) + nt[1] * scalprod(A2, n) + nt[2] * scalprod(A3, n));
      grad[1] = 2 * (np[0] * scalprod(A1, n) + np[1] * scalprod(A2, n) + np[2] * scalprod(A3, n));
    }
    double minPar(int n) const { return (n==0) ? 0    : -PI/2;}
    double maxPar(int n) const { return (n==0) ? 2*PI :  PI/2;}
    P nvec(double theta, double psi) const { // the normal vector, as a function of theta and psi
      return {std::cos(theta) * std::cos(psi),
              std::sin(theta) * std::cos(psi),
              std::sin(psi)};
    }
    P nvec_dtheta(double theta, double psi) const {
      return {-std::sin(theta) * std::cos(psi),
               std::cos(theta) * std::cos(psi),
               0};
    }
      
    P nvec_dpsi(double theta, double psi) const {
      return {-std::cos(theta) * std::sin(psi),
              -std::sin(theta) * std::sin(psi),
               std::cos(psi)};
    }

    const double PI = 3.1415926536;    
    static double scalprod(const P& p1, const P& p2) {
      return (p1[0] * p2[0]) + (p1[1] * p2[1]) + (p1[2] * p2[2]);
    };
  };
  std::array<double, 2> seed{1, 1}; // @@ will this seed always work?
  Go::FunctionMinimizer<TmpFunc> fmin(2, TmpFunc {A1, A2, A3}, &seed[0]);
  Go::minimise_conjugated_gradient(fmin);

  const double theta = fmin.getPar(0);
  const double psi   = fmin.getPar(1);

  P n {std::cos(theta) * std::cos(psi),
       std::sin(theta) * std::cos(psi),
       std::sin(psi)};

  // we still need to compute the last component of the result array
  const double d = -TmpFunc::scalprod(n, mean_pt);
  return {n[0], n[1], n[2], d};
  
}
  
// ----------------------------------------------------------------------------
template<typename P>
std::array<double, 4> fit_loop_to_plane(const P* const points,
                                        const unsigned int num_points)
// ----------------------------------------------------------------------------
{
  std::array<double, 4> result = fit_points_to_plane(points, num_points);

  // computing the "normal" of the loop
  auto cross = [](const P& p1, const P& p2) { return P {p1[1] * p2[2] - p1[2] * p2[1],
                                                        p1[2] * p2[0] - p1[0] * p2[2],
                                                        p1[0] * p2[1] - p1[1] * p2[0]};};
  P ln {0.0, 0.0, 0.0};
  for (uint i = 0; i != (uint)num_points; ++i) {
    auto& P1 = points[i];
    auto& P2 = points[(i+1) % num_points];
    const P tmp = cross(P1, P2);
    ln[0] = ln[0] + tmp[0]; ln[1] = ln[1] + tmp[1]; ln[2] = ln[2] + tmp[2];
  }

  // take scalar product with plane normal, and check if they point the same way
  if ((ln[0] * result[0] + ln[1] * result[1] + ln[2] * result[2]) < 0) 
    for (uint i= 0; i != 4; ++i)
      result[i] = -result[i];
  return result;
}
  
};

#endif
