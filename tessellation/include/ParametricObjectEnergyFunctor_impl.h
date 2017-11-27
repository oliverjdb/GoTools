#include <iostream> // @@ for debugging only
namespace TesselateUtils {

// ----------------------------------------------------------------------------
ClippedGrid<1> inline
ParamCurveEnergyFunctionTraits::compute_cgrid(const ParamObj pobj,
                                              const BoundaryPolytope& boundary,
                                              const double radius)
// ----------------------------------------------------------------------------  
{
  // ClippedGrid not used in the curve case, so we just return a dummy one.
  return ClippedGrid<1>();
}

// ----------------------------------------------------------------------------  
inline bool
close_to_surface_boundary(const ParamSurfaceEnergyFunctionTraits::ParamObj pobj,
                          std::vector<Point3D> bnd3D,
                          const double radius,
                          const double upar,
                          const double vpar)
// ----------------------------------------------------------------------------  
{
  const auto go_pt(pobj->point(upar, vpar));
  const Point3D pt {go_pt[0], go_pt[1], go_pt[2]};
  uint dummy; 
  return (dist2(pt,
                closest_point_on_loop(pt,
                                      &bnd3D[0],
                                      (uint)bnd3D.size(),
                                      dummy)) < radius*radius) ?
    CLOSE_INSIDE: FAR_INSIDE;
}
                       
// ----------------------------------------------------------------------------
inline void classify_inside_surface_cells(ClippedGrid<2>& cgrid,
              const ParamSurfaceEnergyFunctionTraits::ParamObj pobj,
              const ParamSurfaceEnergyFunctionTraits::BoundaryPolytope& boundary,
              const double radius)
// ----------------------------------------------------------------------------
{
  // first, initialize all inside cells to FAR_INSIDE
  std::transform(cgrid.type.begin(), cgrid.type.end(), cgrid.type.begin(),
                 [](ClippedDomainType type) {
                   return (type == CLOSE_INSIDE) ? FAR_INSIDE : type;});

  // Then, compute the actual 3D boundary
  std::vector<Point3D> bnd3D(boundary.num_bpoints);
  transform(boundary.bpoints, boundary.bpoints + boundary.num_bpoints, bnd3D.begin(),
            [&] (const Point2D& uv) {
              const auto go_point(pobj->point(uv[0], uv[1]));
              return Point3D {go_point[0], go_point[1], go_point[2]}; });
  
  // then, compute the status of all interior corner points: whether they are
  // "close" or "far" from the boundary
  enum Status {NOT_SET, CLOSE, FAR};
  const uint cells_x = cgrid.res[0];
  const uint cells_y = cgrid.res[1];

  std::vector<Status> corner_statuses((cells_x + 1) * (cells_y + 1), NOT_SET);
  for (uint j = 0; j != cgrid.res[1]; ++j) 
    for (uint i = 0; i != cgrid.res[0]; ++i) 
      if (cgrid.type[j * cgrid.res[0] + i] == FAR_INSIDE) 
        // ensure that the status of this cell's corner points have been computed
        for (uint jj = j; jj != j+2; ++jj)
          for (uint ii = i; ii != i+2; ++ii)
            if (corner_statuses[jj * (cgrid.res[0]+1) + ii] == NOT_SET) 
              corner_statuses[jj * (cgrid.res[0]+1) + ii] =
                close_to_surface_boundary(pobj, bnd3D, radius,
                                          cgrid.bbox[0] + ii * cgrid.cell_len[0],
                                          cgrid.bbox[2] + jj * cgrid.cell_len[1]) ?
                CLOSE : FAR;
  
  // classify interior cells
  for (uint j = 0; j != cgrid.res[1]; ++j) 
    for (uint i = 0; i != cgrid.res[0]; ++i) 
      if (cgrid.type[j * cgrid.res[0] + i] == FAR_INSIDE) 
        for (uint jj = j; jj != j+2; ++jj)
          for (uint ii = i; ii != i+2; ++ii)
            if (corner_statuses[jj * (cgrid.res[0] + 1) + ii] == CLOSE)
              // if any of the cell's corners is close to the boundary, the cell
              // is defined to be as well
              cgrid.type[j * cgrid.res[0] + i] = CLOSE_INSIDE;
}
  
// ----------------------------------------------------------------------------
ClippedGrid<2> inline
ParamSurfaceEnergyFunctionTraits::compute_cgrid(const ParamObj pobj,
                                                const BoundaryPolytope& boundary,
                                                const double radius)
// ----------------------------------------------------------------------------  
{
  uint num_cells_x = 40; // @@ hardcoded for now
  uint num_cells_y = 40; // @@ hardcoded for now
  auto cgrid = clip_grid_polygon_2D(boundary.bpoints,
                                    boundary.num_bpoints,
                                    radius,
                                    num_cells_x,
                                    num_cells_y);
                                    
  // INTERSECTED and OUTSIDE are correctly classified, but the above routine
  // cannot distinguish between CLOSE_INSIDE and FAR_INSIDE since we work with a
  // parametrized geometry (distances cannot be directly measured in parameter
  // plane).  We therefore have to call the following function to distinguish
  // between CLOSE_INSIDE and FAR_INSIDE
  
  classify_inside_surface_cells(cgrid, pobj, boundary, radius);
  
  //fill(cgrid.type.begin(), cgrid.type.end(), UNDETERMINED);
  return cgrid;
}

// ----------------------------------------------------------------------------
ClippedGrid<3> inline
ParamVolumeEnergyFunctionTraits::compute_cgrid(const ParamObj pobj,
                                               const BoundaryPolytope& boundary,
                                               const double radius)
// ----------------------------------------------------------------------------  
{
  assert(false);
  return ClippedGrid<3>();
}
  
// ----------------------------------------------------------------------------
template<typename Traits> double
ParametricObjectEnergyFunctor<Traits>::operator()(const double* const arg) const 
// ----------------------------------------------------------------------------
{
  if (!use_cached(arg))
    update_cache(arg);
  return cached_result_.val;
}

// ----------------------------------------------------------------------------
template<typename Traits> void
ParametricObjectEnergyFunctor<Traits>::grad(const double* arg, double* grad) const
// ----------------------------------------------------------------------------
{
  if (!use_cached(arg))
    update_cache(arg);
  const double* const dp = (const double* const)&cached_result_.der[0];
  std::copy(dp, dp + Dim * np_, grad);
}

// ----------------------------------------------------------------------------  
template<typename Traits> bool
ParametricObjectEnergyFunctor<Traits>::use_cached(const double* const arg) const
// ----------------------------------------------------------------------------
{
  return (cached_arg_.size() > 0) &&
         (std::equal(arg, arg + Dim * np_, &cached_arg_[0]));
}

// ----------------------------------------------------------------------------
template<typename Traits> double
ParametricObjectEnergyFunctor<Traits>::minPar(int n) const
// ----------------------------------------------------------------------------  
{
  return bbox_[(n % 2) * 2];
}
  
// ----------------------------------------------------------------------------
template<typename Traits> double
ParametricObjectEnergyFunctor<Traits>::maxPar(int n) const
// ----------------------------------------------------------------------------
{
  return bbox_[(n % 2) * 2 + 1];
}

// ----------------------------------------------------------------------------    
template<> void inline
ParamCurveEnergyFunctor::update_cache(const double* const arg) const
// ----------------------------------------------------------------------------  
{
  if (cached_arg_.empty())
    cached_arg_.resize(Dim * np_);

  std::copy(arg, arg + Dim * np_, &cached_arg_[0]);

  cached_result_ = parametric_curve_energy(pobj_,
                                           minPar(0),
                                           maxPar(0),
                                           &cached_arg_[0],
                                           np_,
                                           radius_);
}

// ----------------------------------------------------------------------------    
template<> void inline
ParamSurfaceEnergyFunctor::update_cache(const double* const arg) const
// ----------------------------------------------------------------------------  
{
  if (cached_arg_.empty())
    cached_arg_.resize(Dim * np_);

  std::copy(arg, arg + Dim * np_, &cached_arg_[0]);

  cached_result_ = parametric_surf_energy(pobj_,
                                          boundary_.bpoints,
                                          boundary_.num_bpoints,
                                          (const Point2D* const) &cached_arg_[0],
                                          np_,
                                          radius_,
                                          &cgrid_);
}
  
}; // end namespace TesselateUtils
