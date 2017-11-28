#ifndef _TESSELATE_PARAMETRIC_VOLUME_H
#define _TESSELATE_PARAMETRIC_VOLUME_H

#include <vector>
#include "common_defs.h"
#include "GoTools/geometry/ParamCurve.h"
#include "MeshXD.h"

namespace TessellateUtils {
  
  // tessellate parametric curve, returning a vector of parameter values
  std::vector<double>
  tessellateParametricCurve(const std::shared_ptr<const Go::ParamCurve> pc,
                           const double vdist);

  Mesh2D
  tessellateParametricSurface(const shared_ptr<const Go::ParamSurface> surf,
                             const Point2D* const bpoints,
                             const uint num_bpoints,
                             const double vdist);
  
  Mesh3D
  tessellateParametricVolume(const shared_ptr<const Go::ParamVolume> pvol,
                            const Point3D* const bpoints,
                            const Triangle* const tris,
                            const uint num_tris,
                            const double vdist);
  
                                               
};

#endif
