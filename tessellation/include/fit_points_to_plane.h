#ifndef _FIT_POINTS_TO_PLANE_H
#define _FIT_POINTS_TO_PLANE_H

#include <array>

namespace TesselateUtils {


// Fit a minimum of 3 points to a plane in 3D space.  The returned value represent the four
// coefficients in the plane equation ax + by + cz + d = 0.  The first three coefficients,
// [a,b,c] represent the plane normal, and have been normalized such that a^2 + b^2 + c^2 = 1.
// This means that the distance from any point (x1, y1, z1) to the plane is given by
// abs(a x1 + b y1 + c z1 + d).
template<typename P>
std::array<double, 4> fit_points_to_plane(const P* const points,
                                          const unsigned int num_points);

// As above, but considers the points to lie on a loop, and chooses orientation of normal to
// ensure loop becomes counterclockwise oriented in plane
template<typename P>
std::array<double, 4> fit_loop_to_plane(const P* const points,
                                        const unsigned int num_points);
};

#include "fit_points_to_plane_impl.h"


#endif
