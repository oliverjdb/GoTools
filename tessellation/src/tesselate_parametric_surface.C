#include <cmath>
#include "ParametricObjectEnergyFunctor.h"
#include "GoTools/utils/GeneralFunctionMinimizer.h"
#include "tesselate_parametric_volume.h"
#include "tesselate_utils.h"
#include "triangulate_domain.h"
#include "ParametricObjectEnergyFunctor.h"

using namespace TesselateUtils;
using namespace std;
using namespace Go;

namespace {

typedef shared_ptr<const ParamSurface> SurfPtr;

vector<Point2D> init_startpoints(const SurfPtr sp,
                                 const Point2D* const bpoints,
                                 const uint num_bpoints,
                                 const double vdist);
void optimize_interior_points(const SurfPtr sp,
                              const Point2D* const bpoints,
                              const uint num_bpoints,
                              Point2D* const ipoints,
                              const uint num_ipoints,
                              const double vdist);
  
}; // end anonymous namespace

namespace TesselateUtils {

  // ============================================================================
Mesh2D tesselateParametricSurface(const SurfPtr surf,
                                  const Point2D* const bpts,
                                  const uint num_bpts,
                                  const double vdist)
// ============================================================================
{
  // choosing an initial set of interior points which can later be moved around
  // to optimal locations.  
  vector<Point2D> ipoints = init_startpoints(surf, bpts, num_bpts, vdist);

  //ipoints = vector<Point2D> {{-1.0, 2.0}, {0.5, 0.5}, {0.6, 0.5}, {0.6, 0.55}}; // @@@
  assert(polygon_area(bpts, num_bpts) > 0);
  
  // optimizing position of interior points
  if (!ipoints.empty())
    optimize_interior_points(surf, bpts, num_bpts,
                            &ipoints[0], (uint)ipoints.size(), vdist * 2.0);

  // triangulating points
  vector<Point2D> points(bpts, bpts + num_bpts);
  points.insert(points.end(), ipoints.begin(), ipoints.end());

  vector<Point3D> points3D(points.size());
  transform(points.begin(), points.end(), points3D.begin(), [surf] (const Point2D& uv) {
      Go::Point p;
      surf->point(p, uv[0], uv[1]);
      return Point3D {p[0], p[1], p[2]};
    });

  vector<Point3D> bpoint_normals(num_bpts);
  transform(bpts, bpts + num_bpts, bpoint_normals.begin(), [surf] (const Point2D& uv) {
      Go::Point result;
      surf->normal(result, uv[0], uv[1]);
      return Point3D {result[0], result[1], result[2]};
    });
  
  // const auto tris = triangulate_domain(&points[0], num_bpts,
  //                                      (uint)points.size(), 3*pardist); // 100 * pardist
  const auto tris = triangulate_2D_manifold_in_3D(&points3D[0],
                                                  num_bpts,
                                                  (uint)points3D.size(),
                                                  3*vdist,
                                                  &bpoint_normals[0],
                                                  &points[0]);
  return {points, tris};
  
}

  
}; // end namespace TesselateUtils

namespace {

// ----------------------------------------------------------------------------  
void optimize_interior_points(const SurfPtr sp,
                              const Point2D* const bpoints,
                              const uint num_bpoints,
                              Point2D* const ipoints,
                              const uint num_ipoints,
                              const double vdist)
// ----------------------------------------------------------------------------  
{
  // setting up function to minimize (energy function)
  const auto efun = ParamSurfaceEnergyFunctor(sp, {bpoints, num_bpoints},
                                              vdist, num_ipoints);

  // setting up function minimizer
  Go::FunctionMinimizer<ParamSurfaceEnergyFunctor>
    funcmin(num_ipoints * 2, efun, (double* const)&ipoints[0], 1e-1 * vdist); // @@ tolerance?1e-1

  // do the minimization
  const double STOPTOL = 1e-4; // @@ is this always sufficiently good?  (Default
                               // in 'minimise_conjugated_gradient' is machine
                               // precision) for this tolerance. 1e-4
  Go::minimise_conjugated_gradient(funcmin, STOPTOL);

  // copying results back.  The cast is based on the knowledge that Points2D
  // consists of POD and that points are stored contiguously in memory
  // (as x1,y1,x2,y2, ...)
  double* const target = (double* const) ipoints;
  std::copy(funcmin.getPar(), funcmin.getPar() + funcmin.numPars(), target);
  
}
  
// ----------------------------------------------------------------------------
inline double compute_surf_param_area(SurfPtr sp)
// ----------------------------------------------------------------------------
{
  const RectDomain d = sp->containingDomain();
  return (d.umax() - d.umin()) * (d.vmax() - d.vmin());
}

// ----------------------------------------------------------------------------
inline double estimate_surf_area(SurfPtr sp)
// ----------------------------------------------------------------------------
{
  // @@ As far as I understand the implementation of estimateSfSize, the area
  // estimated here is tied to the (rectangular) 'containingDomain()', and not
  // the actual trimmed surface.  This is good, since we want a surface area
  // associated with the 'containingDomain()'.
  double u_size, v_size;
  sp->estimateSfSize(u_size, v_size);
  return u_size * v_size;
}
  
// ----------------------------------------------------------------------------
vector<Point2D> init_startpoints(const SurfPtr sp,
                                 const Point2D* const bpoints,
                                 const uint num_bpoints,
                                 const double vdist)
// ----------------------------------------------------------------------------
{
  const double poly_area = polygon_area(bpoints, num_bpoints);

  const double surf_param_area = compute_surf_param_area(sp);
  const double surf_area_estimate = estimate_surf_area(sp);
  const double param_surf_ratio = surf_param_area / surf_area_estimate;

  // estimate average distance in parameter plane
  const double pardist = vdist * sqrt(param_surf_ratio);

  // compute approx. number of points needed to cover whole surface
  const uint N = (uint)floor(poly_area / (2 * pardist * pardist));
  if (N==0)
    return vector<Point2D>(); // no point to insert
  
  // inserting the required number of points inside the bounding box of the
  // clipping polygon
  const array<double, 4> bbox = bounding_box_2D(bpoints, num_bpoints);
  const double bbox_lx = bbox[1] - bbox[0];
  const double bbox_ly = bbox[3] - bbox[2];

  const bool lx_smallest = bbox_lx < bbox_ly;
  const double lmin = (lx_smallest) ? bbox_lx : bbox_ly;
  const double lmax = (lx_smallest) ? bbox_ly : bbox_lx;
  const uint n1 = (uint)ceil(sqrt(N * lmin / lmax));
  const uint n2 = (uint)ceil(N / n1);
  const uint nx = (lx_smallest) ? n1 : n2;
  const uint ny = (lx_smallest) ? n2 : n1;

  // constructing regular grid of points within the bounding box
  vector<Point2D> result = generate_grid_2D(Point2D {bbox[0], bbox[2]},
                                            Point2D {bbox[1], bbox[3]},
                                            nx, ny, false);

  // add perturbation to avoid 'symmetry locking' during the optimization stage
  const double PM = 0.3 * pardist;
  transform(result.begin(), result.end(), result.begin(), 
	    [PM] (Point2D p) { return Point2D {p[0] + random_uniform(-PM, PM),
		                               p[1] + random_uniform(-PM, PM)};});
  return result;
}

// // ----------------------------------------------------------------------------
// vector<Point2D> init_startpoints(const SurfPtr sp,
//                                  const Point2D* const bpoints,
//                                  const uint num_bpoints,
//                                  const double vdist,
//                                  double& pardist)
// // ----------------------------------------------------------------------------
// {
//   const double poly_area = polygon_area(bpoints, num_bpoints);

//   // determine area of _unbounded_ surface (we need the area corresponding to
//   // the use of the full parameter domain)
//   const double surf_param_area = compute_surf_param_area(sp);
//   const double surf_area_estimate = estimate_surf_area(sp);

//   // compute approx. number of points needed to cover whole surface
//   const uint N = (uint)ceil(2 * surf_area_estimate / (vdist * vdist)); 

//   // estimating number of points needed within the bounding box of the clipping
//   // polygon
//   const array<double, 4> bbox = bounding_box_2D(bpoints, num_bpoints);
//   const double bbox_lx = bbox[1] - bbox[0];
//   const double bbox_ly = bbox[3] - bbox[2];
//   const uint Nbox = (uint)ceil(N * surf_param_area / poly_area); // approx. in bounding box
//   const bool lx_smallest = bbox_lx < bbox_ly;
//   const double lmin = (lx_smallest) ? bbox_lx : bbox_ly;
//   const double lmax = (lx_smallest) ? bbox_ly : bbox_lx;
//   const uint n1 = (uint)ceil(sqrt(Nbox * lmin / lmax));
//   const uint n2 = (uint)ceil(Nbox / n1);
//   const uint nx = (lx_smallest) ? n1 : n2;
//   const uint ny = (lx_smallest) ? n2 : n1;

//   // constructing regular grid of points within the bounding box
//   const vector<Point2D> gridpoints =
//     generate_grid_2D(Point2D {bbox[0], bbox[2]},
// 		     Point2D {bbox[1], bbox[3]},
// 		     nx, ny);
  
//   // keeping the points that fall within the original polygon
//   const double TOL_FAC = 0.25 * vdist; // do not keep points too close to the boundary
//   vector<Point2D> result = inpolygon(&gridpoints[0], (unsigned int)gridpoints.size(),
// 				     bpoints,  num_bpoints, TOL_FAC);

//   // add perturbation to avoid 'symmetry locking' during the optimization stage
//   const double PM = 1.0e-1 * min(bbox_lx/nx, bbox_ly/ny); // perturbation magnitude
//   transform(result.begin(), result.end(), result.begin(), 
// 	    [PM] (Point2D p) { return Point2D {p[0] + random_uniform(-PM, PM),
// 		                               p[1] + random_uniform(-PM, PM)};});

//   // estimating parametric distance
//   pardist = min(bbox_lx / nx, bbox_ly / ny);
//   return result;
// }

  
};
