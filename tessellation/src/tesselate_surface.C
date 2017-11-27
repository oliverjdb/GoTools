#include <assert.h>
#include "tesselate_surface.h"
#include "tritools.h"
#include "make_spacing_funs.h"
#include "GoTools/utils/GeneralFunctionMinimizer.h"

using namespace Go;
using namespace std;
using namespace TriTools;

namespace {

class FunctorWrapper
{
public:
  FunctorWrapper(const tuple<RnToRFunction, RnToRnFunction>& fun,
		   unsigned int dim,
		   double minpar,
		   double maxpar) :
    fun_(fun), dim_(dim), minpar_(minpar), maxpar_(maxpar)
  {};

  double operator()(const double* arg) const {
    return get<0>(fun_)(arg, (int)dim_);
  }
  void grad(const double* arg, double* grad) const {
    auto val_and_jac = get<1>(fun_)(arg, dim_);
    copy(val_and_jac.value.begin(), val_and_jac.value.end(), grad);
  }
  double minPar(int n) const { return minpar_; }
  double maxPar(int n) const { return maxpar_; }
  
private:
  const tuple<RnToRFunction, RnToRnFunction> fun_;
  const unsigned int dim_;
  const double minpar_;
  const double maxpar_;
};

  

// ----------------------------------------------------------------------------
vector<double> define_parvec(double minparam, double maxparam, unsigned int num_intparams) {
// ----------------------------------------------------------------------------
  vector<double> result(num_intparams+2);
  result.front() = minparam;
  result.back()  = maxparam;
  for (size_t i = 1; i != num_intparams+1; ++i)
    result[i] = (maxparam-minparam)/(num_intparams + 1) * double(i);
  return result;
}

// ----------------------------------------------------------------------------
vector<Point> setup_parpoint_grid(const vector<double>& u, const vector<double>& v)
// ----------------------------------------------------------------------------
{
  vector<Point> result;
  for (auto vpar : v)
    for (auto upar : u)
      result.push_back({upar, vpar});
  return result;
}
  
}; 


namespace Go {

  // ----------------------------------------------------------------------------
SurfaceTriangulation tesselate_surface(const ParamSurface& ps,
				       const unsigned int num_internal_points_u,
				       const unsigned int num_internal_points_v,
				       const double margin_boundary,
				       const vector<vector<Point>>& bounding_curves)
// ----------------------------------------------------------------------------
{
  // Regularly sample surface
  const RectDomain dom = ps.containingDomain();
  const vector<double> u = define_parvec(dom.umin(), dom.umax(), num_internal_points_u);
  const vector<double> v = define_parvec(dom.vmin(), dom.vmax(), num_internal_points_v);

  vector<Point> parpoints = setup_parpoint_grid(u, v);
  
  // Clip surface against eventual bounding curves
  if (bounding_curves.size() > 0) {
    vector<int> keep = points_inside_loops(bounding_curves, parpoints, margin_boundary);
    vector<Point> parpoints_new;
    for (size_t i = 0; i != parpoints.size(); ++i)
      if (keep[i])
	parpoints_new.push_back(parpoints[i]);
    parpoints.swap(parpoints_new);
  }
  
  // Create triangulation
  SurfaceTriangulation tri = triangulate_with_boundaries(parpoints, bounding_curves);

  // optimize point positions
  const double TOL = 2*sqrt(numeric_limits<double>::epsilon()); // lowest accetable tol for
                                                                // the conjugated gradient
                                                                // fallback algorithm.
  const int ITER = 100;
  const auto fun = make_triang_spacing_fun(tri);

  vector<double> init_par(tri.num_interior_nodes*2, 0);
  for (size_t i = 0; i != tri.num_interior_nodes; ++i) {
    init_par[2*i] = tri.uv[i][0];
    init_par[2*i+1] = tri.uv[i][1];
  }

  UpdateFun updater = [&] (double* x, const double* const dx) {

    const int dim = tri.num_interior_nodes * 2;
    const double MAX_STEP = 0.3;
    const double dx_max = max (*max_element(dx, dx + dim),
			       -*min_element(dx, dx+dim));
    const double fac = (dx_max < MAX_STEP) ? 1 : MAX_STEP/dx_max;
    
    for (size_t i = 0; i != (size_t)tri.num_interior_nodes*2; ++i)
      x[i] += fac * dx[i];
  };

  
  FunctorWrapper funwrap(fun, tri.num_interior_nodes * 2, 0, 1);
  FunctionMinimizer<FunctorWrapper> f_minimizer((int)tri.num_interior_nodes * 2, funwrap, &init_par[0], 3e-8);
  minimise_conjugated_gradient(f_minimizer);
  
  for (size_t i = 0; i != (size_t)tri.num_interior_nodes; ++i) {
    tri.uv[i][0] = f_minimizer.getPar()[2*i];
    tri.uv[i][1] = f_minimizer.getPar()[2*i+1];
  }

  tri = triangulate_with_boundaries(tri.uv,  bounding_curves);

  //vector<double> new_par = find_root(get<1>(fun), &init_par[0], init_par.size(), updater, TOL, ITER);

  // for (size_t i = 0; i != (size_t)tri.num_interior_nodes; ++i) {
  //   tri.uv[i][0] = new_par[2*i];
  //   tri.uv[i][1] = new_par[2*i+1];
  // }
  return tri;
}


  
}; // end namespace Go;
