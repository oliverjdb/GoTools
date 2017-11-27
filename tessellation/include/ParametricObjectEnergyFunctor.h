#ifndef _PARAMETRIC_OBJECT_ENERGY_FUNCTOR_H
#define _PARAMETRIC_OBJECT_ENERGY_FUNCTOR_H

#include <vector>
#include <memory>
#include "common_defs.h"
#include "tesselate_utils.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/ParamSurface.h"
#include "parametric_object_energies.h"

namespace TesselateUtils {

struct ParamCurveEnergyFunctionTraits {
  static const int Dim = 1; // a 1-manifold (although in 3D space)
  typedef const shared_ptr<const Go::ParamCurve> ParamObj;
  struct BoundaryPolytope {
    const Point1D* const bpoints; // should point to exactly two 1D points
    static const uint num_bpoints = 2;
  };
  static ClippedGrid<1> compute_cgrid(const ParamObj pobj,
                                      const BoundaryPolytope& boundary,
                                      const double radius);
};

struct ParamSurfaceEnergyFunctionTraits {
  static const int Dim = 2; // a 2-manifold (although in 3D space)
  typedef const shared_ptr<const Go::ParamSurface> ParamObj;
  struct BoundaryPolytope {
    const Point2D* const bpoints;
    const uint num_bpoints;
  } ;
  static ClippedGrid<2> compute_cgrid(const ParamObj pobj,
                                      const BoundaryPolytope& boundary,
                                      const double radius);
  
};

struct ParamVolumeEnergyFunctionTraits {
  static const int Dim = 3; // a 3-manifold in 3D-space
  typedef const shared_ptr<const Go::ParamVolume> ParamObj;
  struct BoundaryPolytope {
    const Point3D* const bpoints;
    const uint num_bpoints;
    const Triangle* const btris;
    const uint num_btris;
  };
  static ClippedGrid<3> compute_cgrid(const ParamObj pobj,
                                      const BoundaryPolytope& boundary,
                                      const double radius);
};
  
  
// ============================================================================
template<typename ParamObjectTraits> class ParametricObjectEnergyFunctor  
// ============================================================================
{
  static const int Dim = ParamObjectTraits::Dim;
  typedef typename ParamObjectTraits::ParamObj ParamObj;
  typedef typename ParamObjectTraits::BoundaryPolytope BoundaryPolytope;
  
public:
  ParametricObjectEnergyFunctor(const ParamObj pobj,
                                const BoundaryPolytope& boundary, 
                                double radius,
                                unsigned int num_points)
    : pobj_(pobj), boundary_(boundary), radius_(radius), np_(num_points),
      bbox_(compute_bounding_box(boundary_.bpoints, boundary.num_bpoints)),
      cgrid_(ParamObjectTraits::compute_cgrid(pobj_, boundary_, radius_)) {}

  double operator() (const double* const arg) const;
  void grad(const double* arg, double* grad) const;
  double minPar(int i) const;
  double maxPar(int i) const;
  
private:
  const ParamObj pobj_;                  // parametric object
  const BoundaryPolytope boundary_;      // clipping boundary 
  const double radius_;                  // radius of energy function
  const uint np_;                        // number of unknown points
  const std::array<double, 2*Dim> bbox_; // bounding box
  const ClippedGrid<Dim> cgrid_;         // precomputed classification of
                                         // subdivided domain parts, according
                                         // to their relationship to the polygon
                                         // boundary (inside, outside, etc.).
                                         // Used to improve computational
                                         // efficiency.

  mutable ValAndDer<PointXD<Dim>> cached_result_;
  mutable std::vector<double> cached_arg_;

  bool use_cached(const double* const arg) const;
  void update_cache(const double* const arg) const;
  
}; // end ParametricObjectEnergyFunctor

typedef
ParametricObjectEnergyFunctor<ParamCurveEnergyFunctionTraits>
ParamCurveEnergyFunctor;

typedef
ParametricObjectEnergyFunctor<ParamSurfaceEnergyFunctionTraits>
ParamSurfaceEnergyFunctor;

typedef
ParametricObjectEnergyFunctor<ParamVolumeEnergyFunctionTraits>
ParamVolumeEnergyFunctor;

} // end namespace TesselateUtils

#include "ParametricObjectEnergyFunctor_impl.h"

#endif
