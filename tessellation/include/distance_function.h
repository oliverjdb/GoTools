#ifndef _DISTANCE_FUNCTION_H
#define _DISTANCE_FUNCTION_H

#include <functional>
#include "GoTools/geometry/ParamSurface.h"
#include "common_defs.h"

namespace TesselateUtils
{

  typedef std::function<double(const Point2D&, const Point2D&, double* grad)> SquaredDistanceFun2D;

  static const SquaredDistanceFun2D default_squared_distance_fun_2D =
    [] (const Point2D& p1, const Point2D& p2, double* grad)
  {
    const double dx = p1[0] - p2[0];
    const double dy = p1[1] - p2[1];
    const double dist2 = dx*dx + dy*dy;

    if (grad) {
      // with respect to components of p1
      grad[0] = 2 * dx;
      grad[1] = 2 * dy;
      grad[2] = -2 * dx;
      grad[3] = -2 * dy;
    }

    return dist2;
  };

  SquaredDistanceFun2D make_squared_distance_function(const shared_ptr<const Go::ParamSurface> s);
  
};// end namespace TesselateUtils

#endif
