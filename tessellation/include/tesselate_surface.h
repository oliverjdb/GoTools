#ifndef _TESSELATE_SURFACE_H
#define _TESSELATE_SURFACE_H

#include <vector>
#include "GoTools/geometry/ParamSurface.h"
#include "tritools.h"

namespace Go {
  
  // ----------------------------------------------------------------------------
  // Tesselate (triangulate) parametric surface, with nodes obtained as regular
  // samples, whose numbers in the u and v parameter directions are specified
  // directly.  If the surface is bounded, the bounding curves (external and
  // internal) are provided in the 'bounding_curves' vector, and represents the
  // boundaries in the _parametric_ domain.  If the surface is not bounded by
  // curves, this vector should be empty.  'margin_boundary' specifies how close
  // to the boundary internal sample points are allowed to be (in the case the
  // surface is bounded by curves).
  TriTools::SurfaceTriangulation
  tesselate_surface(const Go::ParamSurface& ps,
		    const unsigned int num_internal_points_u,
		    const unsigned int num_internal_points_v,
		    const double margin_boundary,
		    const std::vector<std::vector<Go::Point>>& bounding_curves);
  // ----------------------------------------------------------------------------
  
};

#endif
