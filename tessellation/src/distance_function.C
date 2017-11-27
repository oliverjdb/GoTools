#include <vector>
#include "distance_function.h"

using namespace std;

namespace TesselateUtils {

// ----------------------------------------------------------------------------  
SquaredDistanceFun2D
make_squared_distance_function(const shared_ptr<const Go::ParamSurface> surf)
// ----------------------------------------------------------------------------
{
  vector<Go::Point> pvec1(3), pvec2(3); // avoid having to recreate these each time function is called 
  return [surf, pvec1, pvec2] (const Point2D& p1, const Point2D& p2, double* grad) mutable {

    auto pt1 = pvec1[0]; // for notational simplicity
    auto pt2 = pvec2[0];

    if (grad) { // compute both distance and gradient

      surf->point(pvec1, p1[0], p1[1], 1);  // compute point and directional derivatives
      surf->point(pvec2, p2[0], p2[1], 1);  // compute point and directional derivatives

    } else { // compute oniy distance
      
      surf->point(pt1, p1[0], p1[1]);
      surf->point(pt2, p2[0], p2[1]);
    }

    // compute squared distance value (unrolling scalar product)
    
    const double dx = pt2[0] - pt1[0];
    const double dy = pt2[1] - pt1[1];
    const double dz = pt2[2] - pt1[2];
    
    const double dist2 = dx*dx + dy*dy + dz*dz;

    // if gradient requested, compute as well
    if (grad) 
      for (uint i = 0; i != 2; ++i) {
        // partial derivatives wrt. parameters for point p2
        grad[i] = -2 * dx * pvec1[i+1][0] +
                       dy * pvec1[i+1][1] +
                       dz * pvec1[i+1][2];
        
        // partial derivatives wrt. parameters for point p1
        grad[2+i] = 2 * dx * pvec2[1+i][0] + 
                        dy * pvec2[1+i][1] +
                        dz * pvec2[1+i][2];
      }

    // return square distance
    return dist2;
  };
}
  
}; //end namespace TesselateUtils
