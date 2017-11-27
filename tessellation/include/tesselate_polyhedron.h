#ifndef _TESSELATE_POLYHEDRON_H
#define _TESSELATE_POLYHEDRON_H

#include <vector>
#include <array>
#include <cmath>
#include "common_defs.h"
#include "distance_function.h"
#include "tesselate_utils.h"
#include "triangulate_domain.h"
#include "GoTools/geometry/ParamSurface.h"
#include "MeshXD.h"

namespace TesselateUtils {

template<typename PointXD> 
std::vector<PointXD> tesselateSegment(const PointXD& p1,
					const PointXD& p2,
					const double vdist);

// Boundary polygon should be given counterclockwise
Mesh2D tesselatePolygon2D(const Point2D* const polygon,
			  const unsigned int num_corners,
			  const double vdist,
                          const shared_ptr<const Go::ParamSurface> surf = 0);

// All face triangle should be presented with face normal pointing outwards.
Mesh3D tesselatePolyhedron3D(const Point3D* const bpoints,
                             const unsigned int num_bpoints,
                             const Triangle* const btris,
                             const unsigned int num_tris,
                             const double vdist);


// ========================= TEMPLATE IMPLEMENTATIONS =========================

template<typename PointXD> inline
std::vector<PointXD> tesselateSegment(const PointXD& p1,
				      const PointXD& p2,
				      const double vdist)
{
  const double seg_len = dist(p1, p2);
  const unsigned int num_intervals = (unsigned int)std::ceil(seg_len/vdist);

  return interpolate(p1, p2, num_intervals - 1);
}
};

#endif
