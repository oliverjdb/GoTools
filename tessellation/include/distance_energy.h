#ifndef _DISTANCE_ENERGY_H
#define _DISTANCE_ENERGY_H

#include <array>

namespace TesselateUtils {

// ----------------------------------------------------------------------------
// Energy as a function of distance.  Its support remains within the support of
// R.  Returns energy and derivative.
inline std::array<double, 2> distance_energy(double dist, double R)
// ----------------------------------------------------------------------------
{
  const double tmp = std::max(R-dist, double(0));
  const double tmp2 = tmp*tmp;
  const double R4 = R*R*R*R; // normalizing factor so that E(0) = 1
  return {tmp2 * tmp2/R4, -4 * tmp*tmp*tmp/R4}; // energy and derivative
}

  
  
};

#endif
