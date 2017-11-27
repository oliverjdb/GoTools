#include "make_spacing_funs.h"
#include "GoTools/utils/Point.h"

using namespace std;
using namespace Go;

namespace {

// Evaluate a set of points (and optional derivatives) along the ParamCurve,
// including start and end points.
// par - refers to parameters for _interior_ points
vector<vector<Point>>
compute_points(const ParamCurve& c, const double* const par,
	       unsigned int dim, unsigned int ders = 0);

double    eval_fun(const ParamCurve& c, const double* const par, unsigned int dim);
ValAndJac eval_dfun(const ParamCurve& c, const double* const par, unsigned int dim); 
  
}; // end anonymous namespace



namespace Go {

// ============================================================================
tuple<RnToRFunction, RnToRnFunction> make_curve_spacing_fun(const ParamCurve& c)
// ============================================================================
{
  return tuple<RnToRFunction, RnToRnFunction> {
    [&c] (const double* const par, unsigned int dim) {return eval_fun (c, par, dim);},
    [&c] (const double* const par, unsigned int dim) {return eval_dfun(c, par, dim);}
  };
}
  
}; // end namespace Go

namespace {

// ----------------------------------------------------------------------------
vector<vector<Point>>
compute_points(const ParamCurve& c, const double* const par,
	       unsigned int dim, unsigned int ders)
// ----------------------------------------------------------------------------
{
  vector<vector<Point>> points(dim+2, vector<Point>(ders+1));

  // Computing start and end points
  c.point(points.front(), c.startparam(), ders);
  c.point(points.back(), c.endparam(), ders);
  
  // Computing interior points
  for (size_t i = 0; i != dim; ++i)
    c.point(points[i+1], par[i], ders);

  return points;
}
  
// ----------------------------------------------------------------------------
double eval_fun(const ParamCurve& c, const double* const par, unsigned int dim)
// ----------------------------------------------------------------------------
{
  const auto points = compute_points(c, par, dim, 0);

  // computing square distances
  double res = 0;
  for (auto p1 = points.begin(), p2 = p1+1; p2 != points.end(); ++p1, ++p2)
    res = res + (*p1)[0].dist2((*p2)[0]);
  
  return res;
}

// ----------------------------------------------------------------------------
ValAndJac eval_dfun(const ParamCurve& c, const double* const par, unsigned int dim)
// ----------------------------------------------------------------------------
{
  const auto points = compute_points(c, par, dim, 2);
  
  // computing first derivative
  vector<double> d1(dim);
  for (size_t i = 1; i != dim+1; ++i) {
    d1[i-1] = (4 * points[i][0] - 2 * (points[i+1][0] + points[i-1][0])) * points[i][1];
  }
  
  // computing second derivative
  vector<double> d2_diag(dim), d2_subdiag(dim-1);
  for (size_t i = 1; i != dim+1; ++i) {
    d2_diag[i-1] =
      (4 * points[i][0] - 2 * (points[i+1][0] + points[i-1][0])) * points[i][2] +
      4 * points[i][1].length2();
    if (i>1) d2_subdiag[i-2]   = -2 * points[i-1][1] * points[i][1];
  }

  // constructing Jacobian
  vector<double> d2(dim * dim);
  for (size_t i = 0; i != dim; ++i) 
    d2[i * dim + i] = d2_diag[i]; // filling in diagonal
  for (size_t i = 0; i != dim-1; ++i) {
    d2[i*dim + i + 1] = d2_subdiag[i];
    d2[(i+1)*dim + i] = d2_subdiag[i];
  }
  return {d1, d2};
};
  
}; // end anonymous namespace
	    
