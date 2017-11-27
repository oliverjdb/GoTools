#include <iostream> // for debugging
#include <fstream> // for debugging
#include <iterator> // for debugging
#include <assert.h>
#include <cmath>
#include <algorithm>
#include "clip_grid.h"
#include "tesselate_utils.h"
#include "tesselate_polyhedron.h"
#include "polyhedral_energies.h"
#include "GoTools/utils/GeneralFunctionMinimizer.h"
#include "distance_function.h"

#include <string> // remove

using namespace Go;
using namespace std;
using namespace TesselateUtils; 

namespace {


const double PI = 3.14159265358979323;

// Class needed by the "optimize_interior_point" function for the call to
// Go::minimise_conjugated_gradient.
class PolygonEnergyFunctor
{
public:
  PolygonEnergyFunctor(const Point2D* const polygon,
		       const unsigned int num_corners,
		       const unsigned int num_ipoints,
		       const double radius);

  double operator()(const double* const arg) const;
  void grad(const double* arg, double* grad) const;
  double minPar(int n) const;
  double maxPar(int n) const;
private:

  bool use_cached(const double* const arg) const;
  void update_cache(const double* const arg) const;
  const Point2D* const poly_;
  const unsigned int nc_; // num polygon corners
  const unsigned int ni_; // num interior points
  const double r_; // radius
  const std::array<double, 4> bbox_; // bounding box
  const ClippedGrid<2> cgrid_; // precomputed classification of subdivided
                               // domain parts, according to their relationship
                               // to the polygon boundary (inside, outside,
                               // etc.).  Used to improve computational efficiency.

  mutable ValAndDer<Point2D> cached_result_;
  mutable vector<double> cached_arg_;

};

// Class needed by the "optimize_interior_point" function for the call to
// Go::minimise_conjugated_gradient.
class PolyhedronEnergyFunctor
{
public:
  PolyhedronEnergyFunctor(const Point3D* const bpoints,
                          const unsigned int num_bpoints,
                          const Triangle* const btris,
                          const unsigned int num_btris,
                          const Point3D* const ipoints,
                          const unsigned int num_ipoints,
                          const double radius);

  double operator()(const double* const arg) const;
  void grad(const double* arg, double* grad) const;
  double minPar(int n) const;
  double maxPar(int n) const;
private:

  bool use_cached(const double* const arg) const;
  void update_cache(const double* const arg) const;
  const Point3D* const bpoints_;
  const unsigned int nb_; // num boundary points
  const Triangle* const btris_;
  const unsigned int nt_; // num triangles
  const unsigned int ni_; // num interior points
  const double r_; // radius
  const std::array<double, 6> bbox_; // bounding box
  const ClippedGrid<3> cgrid_; // precomputed classification of subdivided domain parts

  mutable ValAndDer<Point3D> cached_result_;
  mutable vector<double> cached_arg_;
};

// ----------------------------------------------------------------------------  
vector<Point2D> init_startpoints(const Point2D* const polygon,
				 const unsigned int num_corners,
				 const double vdist,
                                 const shared_ptr<const ParamSurface> surf);
// choose some initial startpoints as a basis for subsequent optimization  
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
vector<Point3D> init_startpoints(const Point3D* const bpoints,
                                 const unsigned int num_bpoints,
                                 const Triangle* btris,
                                 const unsigned int num_btris,
                                 const double vdist);
// choose some initial startpoints as a basis for subsequent optimization
// ----------------------------------------------------------------------------
  
// ----------------------------------------------------------------------------
void optimize_interior_points(const Point2D* const polygon,
			      const unsigned int num_corners,
			      Point2D* const ipoints,
			      const unsigned int num_ipoints,
			      const double vdist);
// ----------------------------------------------------------------------------  


// ----------------------------------------------------------------------------    
void optimize_interior_points(const Point3D* bpoints,
                              const unsigned int num_bpoints,
                              const Triangle* const btris,
                              const unsigned int num_btris,
                              Point3D* const ipoints,
                              const unsigned int num_ipoints,
                              const double vdist);
// ----------------------------------------------------------------------------    
  
}; // end anonymous namespace

namespace TesselateUtils {

// ============================================================================
Mesh2D tesselatePolygon2D(const Point2D* const polygon,
			  const unsigned int num_corners,
			  const double vdist,
                          const shared_ptr<const ParamSurface> surf)
// ============================================================================
{
  // Choosing a set of interior points which we can later move around to optimal
  // locations 
  vector<Point2D> ipoints = init_startpoints(polygon, num_corners, vdist, surf); 

  // Optimizing position of interior points
  if (!ipoints.empty())
    optimize_interior_points(polygon, num_corners,
                             &ipoints[0], (uint)ipoints.size(),
                             vdist * 1.5); //1); // * 1.5 @@
  // Triangulating points
  vector<Point2D> points(polygon, polygon + num_corners);
  points.insert(points.end(), ipoints.begin(), ipoints.end());

  const auto tris = triangulate_domain(&points[0],
				       num_corners,
				       (uint)points.size(),
				       3*vdist);
  return {points, tris};
}


  
// ============================================================================
Mesh3D tesselatePolyhedron3D(const Point3D* const bpoints,
                             const unsigned int num_bpoints,
                             const Triangle* const btris,
                             const unsigned int num_btris,
                             const double vdist)
// ============================================================================
{
  // choosing a set of interior points
  vector<Point3D> ipoints = init_startpoints(bpoints, num_bpoints, btris,
                                             num_btris, vdist);
  
  // optimizing position of interior points.  We use a value of 'vdist' slightly
  // higher than what the interpoint distance goal is, to avoid points becoming
  // 'completely disconnected' from each other.
  cout << "Optimizing interior points..." << endl;
  if (!ipoints.empty())
    optimize_interior_points(bpoints, num_bpoints, btris, num_btris,
                             &ipoints[0], (uint)ipoints.size(), vdist*2);

  cout << "Finished optimization of " << ipoints.size() << " internal points." << endl;
  
  // In case boundary has much larger faces than 'vdist', we need to allow for
  // larger distances when constructing the tetrahedrons.  We therefore run a
  // check here, and adjust the 'vdist' parameter upwards if the boundary faces
  // require it.
  vector<double> lengths(num_btris, 0);
  transform(btris, btris + num_btris, lengths.begin(), [bpoints] (const Triangle& t) {
      return max({dist(bpoints[t[0]], bpoints[t[1]]),
                  dist(bpoints[t[1]], bpoints[t[2]]),
                  dist(bpoints[t[2]], bpoints[t[0]])});});
  const double L = *max_element(lengths.begin(), lengths.end());
  
  // constructing tets from points
  cout << "Constructing tets" << endl;
  vector<Point3D> points(bpoints, bpoints + num_bpoints);
  points.insert(points.end(), ipoints.begin(), ipoints.end());
  const auto tets = construct_tets(&points[0], (uint)points.size(),
                                   btris, num_btris, 1.5 * max(vdist,L));
  
  cout << "Finished constructing tets" << endl;
  return{points, tets};
}
};

namespace {

// ----------------------------------------------------------------------------
vector<Point3D> init_startpoints(const Point3D* const bpoints,
                                 const unsigned int num_bpoints,
                                 const Triangle* btris,
                                 const unsigned int num_btris,
                                 const double vdist)
// choose some initial startpoints as a basis for subsequent optimization
// ----------------------------------------------------------------------------
{
  const array<double, 6> bbox = bounding_box_3D(bpoints, num_bpoints);
  const double bbox_lx = bbox[1] - bbox[0];
  const double bbox_ly = bbox[3] - bbox[2];
  const double bbox_lz = bbox[5] - bbox[4];

  // computing amount of points needed to approximately fill the bounding box,
  // where points are approximately equidistant by 'vdist'.
  const int nx = 1.1 * max((int)floor(bbox_lx/vdist), 0);
  const int ny = 1.1 * max((int)floor(bbox_ly/vdist), 0);
  const int nz = 1.1 * max((int)floor(bbox_lz/vdist), 0);

  // constructing regular grid of points within the bounding box
  const vector<Point3D> gridpoints =
    generate_grid_3D(Point3D {bbox[0], bbox[2], bbox[4]},
                     Point3D {bbox[1], bbox[3], bbox[5]},
                     nx, ny, nz);

  // keeping the points that fall within the shell
  const double TOL_FAC = 0.25 * vdist; // do not keep points too close to boundary
  vector<Point3D> result = inside_shell(&gridpoints[0],
                                        (unsigned int) gridpoints.size(),
                                        bpoints,
                                        btris,
                                        num_btris,
                                        TOL_FAC);

  // add perturbation to avoid 'symmetry locking' during the optimization stage
  const double PM = 1.0e-1 * min({bbox_lx/nx, bbox_ly/ny, bbox_lz/nz});
  transform(result.begin(), result.end(), result.begin(), [PM] (const Point3D& p) {
      return Point3D {p[0] + random_uniform(-PM, PM),
                      p[1] + random_uniform(-PM, PM),
                      p[2] + random_uniform(-PM, PM)};});
  return result;
}
   
// ----------------------------------------------------------------------------  
vector<Point2D> init_startpoints(const Point2D* const polygon,
				 const unsigned int num_corners,
				 const double vdist,
                                 const shared_ptr<const ParamSurface> surf)
// choose some initial startpoints as a basis for subsequent optimization  
// ----------------------------------------------------------------------------
{
  const double poly_area = polygon_area(polygon, num_corners);

  // If points are to be optimized across a paramtertic surface, use area of
  // surface to estimate number of points.  Otherwise, points are to be
  // optimized directly across the parameter domain.  In that case, use the area
  // of the polygon itself.
  const double SURF_AREA_TOL = 1e-2;
  const double surf_area = surf ? surf->area(SURF_AREA_TOL) : poly_area;
  const uint N = (uint)ceil(2 * surf_area / (vdist * vdist)); // approximate number of points needed
  
  const array<double, 4> bbox = bounding_box_2D(polygon, num_corners);
  const double bbox_lx = bbox[1] - bbox[0];
  const double bbox_ly = bbox[3] - bbox[2];
  const double bbox_area =  bbox_lx * bbox_ly;

  const uint Nbox = (uint)ceil(N * bbox_area / poly_area); // approx. number of points in bounding box
  const bool lx_smallest = bbox_lx < bbox_ly;
  const double lmin = (lx_smallest) ? bbox_lx : bbox_ly;
  const double lmax = (lx_smallest) ? bbox_ly : bbox_lx;
  const uint n1 = (uint)ceil(sqrt(Nbox * lmin / lmax));
  const uint n2 = (uint)ceil(Nbox / n1);
  const uint nx = (lx_smallest) ? n1 : n2;
  const uint ny = (lx_smallest) ? n2 : n1;

  // computing amount of points needed to approximately cover the bounding box, where points
  // are approximately equidistant by 'vdist'
  // const int nx = max((int)floor(bbox_lx/vdist), 0);
  // const int ny = max((int)floor(bbox_ly/vdist), 0);
  
  // constructing regular grid of points within the bounding box
  const vector<Point2D> gridpoints =
    generate_grid_2D(Point2D {bbox[0], bbox[2]},
		     Point2D {bbox[1], bbox[3]},
		     nx, ny);
  
  // keeping the points that fall within the original polygon
  const double TOL_FAC = 0.25 * vdist; // do not keep points too close to the boundary
  vector<Point2D> result = inpolygon(&gridpoints[0], (unsigned int)gridpoints.size(),
				     polygon,  num_corners, TOL_FAC);

  // add perturbation to avoid 'symmetry locking' during the optimization stage
  const double PM = 1.0e-1 * min(bbox_lx/nx, bbox_ly/ny); // perturbation magnitude
  transform(result.begin(), result.end(), result.begin(), 
	    [PM] (Point2D p) { return Point2D {p[0] + random_uniform(-PM, PM),
		                               p[1] + random_uniform(-PM, PM)};});

  return result;
}

// ----------------------------------------------------------------------------
void optimize_interior_points(const Point3D* bpoints,
                              const unsigned int num_bpoints,
                              const Triangle* const btris,
                              const unsigned int num_btris,
                              Point3D* const ipoints,
                              const unsigned int num_ipoints,
                              const double vdist)
// ----------------------------------------------------------------------------
{
  // setting up function to minimize (energy function):
  auto efun = PolyhedronEnergyFunctor(bpoints,num_bpoints, btris, num_btris,
                                      ipoints, num_ipoints, vdist);

  // setting up function minimizer
  Go::FunctionMinimizer<PolyhedronEnergyFunctor>
    funcmin(num_ipoints * 3, efun, (double* const)&ipoints[0], 1e-1); // @@ tolerance?

  const double STOPTOL = 1e-4;

  // do the minimization
  Go::minimise_conjugated_gradient(funcmin, STOPTOL);

  // copying results back.  The cast is based on the knowledge that Points2D
  // consist of POD and that point coordinates are stored contiguously in memory
  // as (x1, y1, z1, x2, y2, z2, ...).
  double* const target = (double* const) ipoints;
  std::copy(funcmin.getPar(), funcmin.getPar() + funcmin.numPars(), target);
}
  
// ----------------------------------------------------------------------------
void optimize_interior_points(const Point2D* const polygon,
			      const unsigned int num_corners,
			      Point2D* const ipoints,
			      const unsigned int num_ipoints,
			      const double vdist)
// ----------------------------------------------------------------------------  
{
  // Setting up function to minimize (energy function);
  auto efun = PolygonEnergyFunctor(polygon, num_corners, num_ipoints, vdist);
  
  // Setting up function minimizer
  Go::FunctionMinimizer<PolygonEnergyFunctor>
    funcmin(num_ipoints * 2, efun, (double* const)&ipoints[0], 1e-1);//1e-1); // @@ TOLERANCE?

  // do the minimization 
  const double STOPTOL = 1e-4; // @@ is this always sufficiently good?  (Default
                               // in 'minimise_conjugated_gradient' is machine
                               // precision) for this tolerance.
  Go::minimise_conjugated_gradient(funcmin, STOPTOL);

  // copying results back.  The cast is based on the knowledge that Points2D
  // consists of POD and that points are stored contiguously in memory
  // (as x1,y1,x2,y2, ...)
  double* const target = (double* const) ipoints;
  std::copy(funcmin.getPar(), funcmin.getPar() + funcmin.numPars(), target);
}

// ----------------------------------------------------------------------------
PolyhedronEnergyFunctor::PolyhedronEnergyFunctor(const Point3D* const bpoints,
                                                 const unsigned int num_bpoints,
                                                 const Triangle* const btris,
                                                 const unsigned int num_btris,
                                                 const Point3D* const ipoints,
                                                 const unsigned int num_ipoints,
                                                 const double radius)
// ----------------------------------------------------------------------------
  : bpoints_(bpoints), nb_(num_bpoints), btris_(btris), nt_(num_btris),
    ni_(num_ipoints), r_(radius), bbox_(bounding_box_3D(bpoints, num_bpoints)),
    cgrid_(clip_grid_shell_3D(bpoints, num_bpoints, btris, num_btris, radius, 20, 20, 20))
    // @@ 50 hard-coded here
{}

// ----------------------------------------------------------------------------
double PolyhedronEnergyFunctor::operator()(const double* const arg) const
// ----------------------------------------------------------------------------
{
  if (!use_cached(arg)) 
    update_cache(arg);
  return cached_result_.val;
}

// ----------------------------------------------------------------------------
void PolyhedronEnergyFunctor::grad(const double* arg, double* grad) const
// ----------------------------------------------------------------------------
{
  if (!use_cached(arg))
    update_cache(arg);
  const double* const dp = (const double* const)&cached_result_.der[0];
  copy(dp, dp + 3 * ni_, grad);
}

// ----------------------------------------------------------------------------
double PolyhedronEnergyFunctor::minPar(int n) const
// ----------------------------------------------------------------------------
{
  return bbox_[(n % 3) * 2];
}

// ----------------------------------------------------------------------------
double PolyhedronEnergyFunctor::maxPar(int n) const
// ----------------------------------------------------------------------------
{
  return bbox_[(n % 3) * 2 + 1];
}

// ----------------------------------------------------------------------------
bool PolyhedronEnergyFunctor::use_cached(const double* const arg) const
// ----------------------------------------------------------------------------
{
  return (cached_arg_.size() > 0) &&
         (std::equal(arg, arg + 3 * ni_, &cached_arg_[0]));
}

// ----------------------------------------------------------------------------  
void PolyhedronEnergyFunctor::update_cache(const double* const arg) const
// ----------------------------------------------------------------------------
{
  if (cached_arg_.empty())
    cached_arg_.resize(3 * ni_);
  copy(arg, arg + 3 * ni_, &cached_arg_[0]);
  cached_result_ = polyhedron_energy(bpoints_, nb_, btris_, nt_,
                                     (const Point3D* const) &cached_arg_[0],
                                     ni_, r_, &cgrid_);
}

// ----------------------------------------------------------------------------    
PolygonEnergyFunctor::PolygonEnergyFunctor(const Point2D* const polygon,
					   const unsigned int num_corners,
					   const unsigned int num_ipoints,
					   const double radius)
// ----------------------------------------------------------------------------    
  : poly_(polygon), nc_(num_corners), ni_(num_ipoints), r_(radius),
    bbox_(bounding_box_2D(polygon, num_corners)),
    cgrid_(clip_grid_polygon_2D(polygon, num_corners, radius, 80, 80)) // @@ 80 hard-coded here
{
}

// ----------------------------------------------------------------------------    
void PolygonEnergyFunctor::update_cache(const double* const arg) const
// ----------------------------------------------------------------------------    
{
  if (cached_arg_.empty())
    cached_arg_.resize(2*ni_);

  copy(arg, arg + 2 * ni_, &cached_arg_[0]);
  //cout << "updating cache." << endl;
  cached_result_ = polygon_energy(poly_, nc_,
				  (const Point2D* const) &cached_arg_[0],
				  ni_, r_, &cgrid_);
}

// ----------------------------------------------------------------------------    
double PolygonEnergyFunctor::operator()(const double* const arg) const
// ----------------------------------------------------------------------------    
{
  // Very inefficient to have separate methods for function value and function
  // gradient.  We have therefore implemented a caching system here.
  if (!use_cached(arg)) {
    update_cache(arg);
    //cout << "updating fval" << endl;
  } else {
    //cout << "using cached fval" << endl;
  }
  return cached_result_.val;
}

// ----------------------------------------------------------------------------    
void PolygonEnergyFunctor::grad(const double* arg, double* grad) const
// ----------------------------------------------------------------------------
{
  // Very inefficient to have separate methods for function value and function
  // gradient.  We have therefore implemented a caching system here.
  if (!use_cached(arg)) {
    update_cache(arg);
    //cout << "updating derivative" << endl;
  } else {
    //cout << "using cached derivative" << endl;
  }
  const double* const dp = (const double* const)&cached_result_.der[0];
  copy(dp, dp + 2 * ni_, grad);
}

// ----------------------------------------------------------------------------    
double PolygonEnergyFunctor::minPar(int n) const
// ----------------------------------------------------------------------------
{
  return bbox_[(n % 2) * 2];
}

// ----------------------------------------------------------------------------
double PolygonEnergyFunctor::maxPar(int n) const 
// ----------------------------------------------------------------------------  
{
  return bbox_[(n % 2) * 2 + 1];
}

// ----------------------------------------------------------------------------  
bool PolygonEnergyFunctor::use_cached(const double* const arg) const
// ----------------------------------------------------------------------------
{
  return (cached_arg_.size() > 0) && (std::equal(arg, arg + 2 * ni_, &cached_arg_[0]));
}

}; // end anonymous namespace
