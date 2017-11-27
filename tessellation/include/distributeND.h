#ifndef _DISTRIBUTE_ND_H
#define _DISTRIBUTE_ND_H

#include <vector>
#include <functional>

namespace Go
{
  // Distance function signature.  First two arguments are the two N-dimensional points whose
  // distance is to be measured.  Third argument, if nonzero, provides derivative information.
  // Fourth argument, if nonzero, provides second derivative information. Function value is
  // returned.
  typedef std::function<double(const double*, const double*, double*, double*)> DistFun;


  std::vector<double> distribute_points_1D(const double* const pts,
					   const unsigned int num_points,
					   const std::pair<double, double>& bounds,
					   const DistFun& dfun);

  std::vector<double> distribute_points_2D(const double* const pts, // u1, v1, u2, v2, ...
					   const unsigned int num_points,
					   const std::vector<std::vector<double>>& bounds,
					   const DistFun& dfun);
};


#endif

