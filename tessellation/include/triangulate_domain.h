#ifndef _TRIANGULATE_DOMAIN_H
#define _TRIANGULATE_DOMAIN_H

#include <array>
#include <vector>
#include "common_defs.h"

namespace TesselateUtils {

  // Compute a generalized Delaunay triangulation of a not necessarily convex
  // domain with a set of given interior points.  The boundary points and
  // interior points are provided by 'points'.  The boundary points are listed
  // first, and there are 'num_bpoints' of them.  The total number of points
  // (boundary points and interior points) is given by 'tot_num_points'.  It is
  // assumed that all provided interior points lie within the boundary, and that
  // the spacing of the points is such that no triangle edge needs to be longer
  // than 'vdist'.  (If you are unsure, you can set 'vdist' as large as you
  // want, but runtime performance is improved by keeping 'vdist' low.  Boundary
  // points should be provided in clockwise order.  Each triangle is returned as
  // a set of three indices in to the list of points.  These indices are
  // provided in counterclockwise order for each triangle.
  std::vector<Triangle> triangulate_domain(const Point2D* const points,
					   const uint num_bpoints, 
					   const uint tot_num_points,
					   const double vdist);

  // Compute a generalized Delaunay triangulation of a not necessarily convex
  // domain with a set of given interior points.  The points are located in 3D
  // space, but assumed to lie on a 2D manifold.  The boundary points and
  // interior points are provided by 'points'.  The boundary points are listed
  // first, and there are 'num_bpoints' of them.  The total number of points
  // (boundary points and interior points) is given by 'tot_num_points'.  It is
  // assumed that all provided interior points lie within the boundary, and that
  // the spacing of the points is such that no triangle edge needs to be longer
  // than 'vdist'.  (If you are unsure, you can set 'vdist' as large as you
  // want, but runtime performance is improved by keeping 'vdist' low.  Boundary
  // points should be provided in clockwise order.  Each triangle is returned as
  // a set of three indices in to the list of points.  These indices are
  // provided in counterclockwise order for each triangle.  Optionally, the user
  // can provide a pointer to an array of normal vectors to the boundary points
  // - in cases where the surface deviates significantly from that of a plane
  // (e.g. U-shaped), this vector of normals might be indispensible.  The user
  // can also optionally provide a parametrization ('param').  If not provided,
  // the function will compute one itself, but in some extreme cases, this has
  // the potential to introduce self-intersections in the final triangulation,
  // and should therefore be avoided on strongly bent surfaces.
  std::vector<Triangle> triangulate_2D_manifold_in_3D(
                                            const Point3D* const points,
                                            const uint num_bpoints,
                                            const uint tot_num_points,
                                            const double vdist,
                                            const Point3D* const bpoint_normals = 0,
                                            const Point2D* const param = 0); 
  
  // Compute a 3D simplex (Tet) mesh from a closed boundary shell defined by a
  // number of triangles defining the boundary, and a number of points defining
  // the boundary points and internal points.  The boundary triangles are
  // expected to reference the N first points of the 'points' vector, where 'N'
  // is the total number of points on the boundary (in other words, boudnary
  // points should be listed first in the array containing points).  It is
  // assumed that the remaining points lie within the interior of the boundary
  // shell, and that the spacing of the points (both interior and on the shell)
  // is such that no tet have any edge longer than 'vdist'.  ('vdist' can be as
  // large as you want, but runtime performance is improved by keeping 'vdist'
  // low).  All boundary triangles should be provied with outward-facing
  // normals, meaning that the three corners making up the triangle are provided
  // in counterclockwise order.
  std::vector<Tet> construct_tets(const Point3D* const points,
                                  const uint tot_num_points,
                                  const Triangle* btris,
                                  const uint num_btris,
                                  const double vdist);
    
};

#endif
