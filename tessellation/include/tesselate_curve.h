#ifndef _TESSELATE_CURVE_H
#define _TESSELATE_CURVE_H

#include <vector>

#include "GoTools/geometry/ParamCurve.h"

namespace Go {

  // Tesselate curve with a prescribed number of internal points.  A vector of parameter values
  // is returned.
  std::vector<double> tesselate_curve(const Go::ParamCurve& pc,
				      const unsigned int num_internal_points);

  // Tesselate curve with a number of internal points adapted to approximately obtain a
  // prescribed spacing between samples.
  std::vector<double> tesselate_curve(const Go::ParamCurve& pc,
				      const double target_spacing);

  // Estimate the length of a parametric curve, based on sampling.
  double estimate_curve_length(const Go::ParamCurve& pc,
			       const unsigned int num_samples = 100);
  
};

#endif
