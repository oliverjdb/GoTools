#ifndef _MAKE_SPACING_FUN_H
#define _MAKE_SPACING_FUN_H

#include <utility>
#include <functional>

#include "GoTools/geometry/ParamCurve.h"
#include "find_root.h"
#include "tritools.h"

using namespace TriTools;

namespace {
  typedef std::function<double(const double* const, int)> RnToRFunction;
};


namespace Go {
// ----------------------------------------------------------------------------
std::tuple<RnToRFunction, RnToRnFunction>
make_curve_spacing_fun(const Go::ParamCurve& c);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------  
std::tuple<RnToRFunction, RnToRnFunction>
make_triang_spacing_fun(const TriTools::SurfaceTriangulation& tri);
// ----------------------------------------------------------------------------
};

#endif // _MAKE_SPACING_FUN_H
