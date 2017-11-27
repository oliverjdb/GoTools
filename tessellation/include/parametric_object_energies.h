#ifndef _PARAMETRIC_OBJECT_ENERGIES_H
#define _PARAMETRIC_OBJECT_ENERGIES_H

#include <vector>
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "clip_grid.h"
#include "common_defs.h"

namespace TesselateUtils
{

// ============================================================================
ValAndDer<Point1D> parametric_curve_energy(const shared_ptr<const Go::ParamCurve> curve,
                                           const double startpar,
                                           const double endpar,
                                           const double* const ipar,
                                           const uint num_ipar,
                                           const double vdist);
// ============================================================================

// ============================================================================
ValAndDer<Point2D> parametric_surf_energy(const shared_ptr<const Go::ParamSurface> surf,
                                          const Point2D* const poly,
                                          const uint num_corners,
                                          const Point2D* const ipar,
                                          const uint num_ipar,
                                          const double vdist,
                                          const ClippedGrid<2>* const cgrid);
// ============================================================================  
  
}; // end namespace TesselateUtils

#endif
