#include "make_spacing_funs.h"

using namespace std;
using namespace Go;

namespace {
double    eval_fun (const SurfaceTriangulation& tri, const double* const par);
ValAndJac eval_dfun(const SurfaceTriangulation& tri, const double* const par); 

inline double cross_value (const double* const p1, const double* const p2)
{
  return p1[0]*p2[1] - p2[0] * p1[1];
}
  
inline double triangle_area(const double* const p1,
			    const double* const p2,
			    const double* const p3)
{
  return 0.5 * (cross_value(p1, p2) + cross_value(p2, p3) + cross_value(p3, p1));
}

inline double squared(double val) {return val*val;}
  
}; // end anonymous namespace

namespace Go {

// ============================================================================
tuple<RnToRFunction, RnToRnFunction>
make_triang_spacing_fun(const SurfaceTriangulation& tri)
// ============================================================================
{
  return tuple<RnToRFunction, RnToRnFunction> {
    // second parameter just a dummy (remains unused) in the following two functions
    [&tri] (const double* const par, unsigned int dummy) {return eval_fun (tri, par);},
    [&tri] (const double* const par, unsigned int dummy) {return eval_dfun(tri, par);}
  };
}
  
}; // end namespace Go

namespace {

// ----------------------------------------------------------------------------
function<const double*(int)> tri_indexer(const SurfaceTriangulation& tri,
				   const double* const par)
// ----------------------------------------------------------------------------
{
  return
    [&tri, par] (int ix) {return ((ix < (int)tri.num_interior_nodes) ?
                          &par[2*ix] :
			  (tri.uv[ix]).begin());};
}

// ----------------------------------------------------------------------------
  function<void(int, int, bool, bool, double)> jac_accum(size_t num_points,
						      double* const data)
// ----------------------------------------------------------------------------  
{
  return [num_points, data] (int i, int j, bool i_y, bool j_y, double val) {
    if ((size_t) max(i, j) < num_points) {
      const size_t i_ix = 2 * i + int(i_y);
      const size_t j_ix = 2 * j + int(j_y);

      data[i_ix + j_ix * 2*num_points] += val;
    }
  };
}
  
// ----------------------------------------------------------------------------
double eval_fun (const SurfaceTriangulation& tri, const double* const par)
// ----------------------------------------------------------------------------  
{
  double result = 0;
  function<const double*(int)> get_pt = tri_indexer(tri, par);

  for (const auto t : tri.triangles) {
    result += squared(triangle_area(get_pt(t[0]), get_pt(t[1]), get_pt(t[2])));
  }
  
  return result;
}

// ----------------------------------------------------------------------------  
ValAndJac eval_dfun(const SurfaceTriangulation& tri, const double* const par)
// ----------------------------------------------------------------------------  
{
  vector<double> value(tri.num_interior_nodes*2, 0l);
  vector<double> jacobian(value.size() * value.size(), 0l);

  function<const double*(int)> get_pt = tri_indexer(tri, par);
  auto accum = jac_accum(tri.num_interior_nodes, &jacobian[0]);
  
  for (const auto t : tri.triangles) {
    const double A = triangle_area(get_pt(t[0]), get_pt(t[1]), get_pt(t[2]));
    
    for (int i = 0; i != 3; ++i) {
      if (t[i] < (int)tri.num_interior_nodes) {
	
	// assembling GRADIENT
	// d/dx
	value[2 * t[i]]   += A * (get_pt(t[(i+1)%3])[1] - get_pt(t[(i+2)%3])[1]);
	// d/dy
	value[2 * t[i]+1] += A * (get_pt(t[(i+2)%3])[0] - get_pt(t[(i+1)%3])[0]);
      }

      // assembling JACOBIAN
      for (int j = 0; j != 3; ++j) {
	// (d/dx) (d/dx)
	accum(t[i], t[j], false, false, 0.5 *
	      (get_pt(t[(j+1)%3])[1] - get_pt(t[(j+2)%3])[1]) *
	      (get_pt(t[(i+1)%3])[1] - get_pt(t[(i+2)%3])[1]));

	// (d/dy) (d/dy)
	accum(t[i], t[j], true, true, 0.5 *
	      (get_pt(t[(j+2)%3])[0] - get_pt(t[(j+1)%3])[0]) *
	      (get_pt(t[(i+2)%3])[0] - get_pt(t[(i+1)%3])[0]));

	// (d/dy) (d/dx)
 	accum(t[i], t[j], false, true, 0.5 *
 	      (get_pt(t[(j+2)%3])[0] - get_pt(t[(j+1)%3])[0]) *
	      (get_pt(t[(i+1)%3])[1] - get_pt(t[(i+2)%3])[1]) +
	      A * (j==((i+1)%3)) - A * (j==((i+2)%3)));

	// // (d/dx) (d/dy)
	accum(t[i], t[j], true, false, 0.5 *
	      (get_pt(t[(j+1)%3])[1] - get_pt(t[(j+2)%3])[1]) *
	      (get_pt(t[(i+2)%3])[0] - get_pt(t[(i+1)%3])[0]) +
	      A * (j==((i+2)%3)) - A * (j==((i+1)%3)));
      }
    }
  }
  return {value, jacobian};
}
  
}; // end anonymous namespace


