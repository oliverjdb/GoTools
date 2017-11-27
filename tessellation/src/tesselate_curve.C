#include "tesselate_curve.h"
#include "make_spacing_funs.h"
#include "find_root.h"
#include "GoTools/utils/GeneralFunctionMinimizer.h"

#include <fstream>
#include <iterator>
#include <ostream>
#include <limits>
#include <cmath>
#include "GoTools/geometry/PointCloud.h"
#include "tesselate_debug.h"

using namespace std;
using namespace Go;

namespace {
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
  vector<Point> evaluate_points(const ParamCurve& pc, const vector<double>& t)
  // ----------------------------------------------------------------------------
  {
    vector<Point> result(t.size());
    for (size_t i = 0; i != t.size(); ++i)
      pc.point(result[i], t[i]);
    return result;
  }

  // ----------------------------------------------------------------------------
  vector<double> compute_distances(const vector<Point> points)
  // ----------------------------------------------------------------------------
  {
    vector<double> res(points.size()-1);
    for (size_t i = 0; i != points.size()-1; ++i) 
      res[i] = points[i].dist(points[i+1]);
    return res;    
  }
  // ----------------------------------------------------------------------------  
  vector<double> adjust_parameters_to_geometry(const ParamCurve& pc, const vector<double>& t)
  // ----------------------------------------------------------------------------
  {

    vector<double> res(t); 

    // evaluate all points
    const vector<Point> points = evaluate_points(pc, t);
    const vector<double> distances = compute_distances(points);

    // If points were equidistant, we estimate the constant distance to be the mean of current
    // distances
    const double edist = accumulate(distances.begin(), distances.end(),
				    0.0)/(double)distances.size();

    // recompute internal parameter values based on point distances
    size_t int_ix = 0; // interval in which we search for current point
    double cur_pos = 0; // current position
    for (size_t par_ix = 1; par_ix != t.size()-1; ++par_ix) {
      const double cur_dist = (double)par_ix * edist; // estimated distance of current point from start point

      // in which interval do we find cur_dist?
      while (cur_pos + distances[int_ix] < cur_dist) 
	cur_pos += distances[int_ix++];

      // now, our target distance should be found between cur_pos and cur_pos + distances[int_ix]
      const double alpha = (cur_dist - cur_pos) / distances[int_ix];

      // setting adjusted parameter value
      res[par_ix] = t[int_ix] + alpha * (t[int_ix+1] - t[int_ix]);
    }
    return res;
  }

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

  
}; // end anonymous namespace

namespace Go {

// ----------------------------------------------------------------------------  
double estimate_curve_length(const ParamCurve& pc, const unsigned int num_samples)
// ----------------------------------------------------------------------------
{
  const vector<double> par =
    adjust_parameters_to_geometry(pc, define_parvec(pc.startparam(),
						    pc.endparam(),
						    num_samples-2));
  Point p, pnew;
  double accum = 0;
  pc.point(p, par[0]);
  for (size_t i = 1; i != par.size(); ++i) {
    pc.point(pnew, par[i]);
    accum += pnew.dist(p);
    p.swap(pnew);
  }
  return accum;
}

// ----------------------------------------------------------------------------  
vector<double> tesselate_curve(const ParamCurve& pc,
			       const double target_spacing)
// ----------------------------------------------------------------------------  
{
  const double curve_len = estimate_curve_length(pc);
  const unsigned int internal_pts = (unsigned int) (curve_len / target_spacing) - 1;
  cout << "internal points: " << internal_pts << endl;
  
  return (internal_pts > 0) ? tesselate_curve(pc, internal_pts) : vector<double>();
}
  
// ----------------------------------------------------------------------------  
vector<double> tesselate_curve(const ParamCurve& pc, const unsigned int num_internal_points)
// ----------------------------------------------------------------------------  
{
  const double TOL = 2*sqrt(numeric_limits<double>::epsilon()); // lowest accetable tol for
                                                              // the conjugated gradient
                                                              // fallback algorithm.
  const int ITER = 100;
  // define a parameter vector with regularly spaced parameters in the domain
  vector<double> p = define_parvec(pc.startparam(), pc.endparam(), num_internal_points);

  //store_points_and_curve(pc, &p[0], (unsigned int)p.size(), "before.g2"); 
			        
  // adjust regularly spaced parameters to better reflect geometry of curve (and
  // thereby hopefully get a better initial guess)
  p = adjust_parameters_to_geometry(pc, p);

  //store_points_and_curve(pc, &p[0], (unsigned int)p.size(), "after.g2");
  
  // define spacing function, and use (modified) Newton-Rhapson to search for a
  // zero of the gradient
  const auto fun = make_curve_spacing_fun(pc);

  // ------------------------- specify updater function -------------------------

  UpdateFun updater = [&] (double* x, const double* const dx) {

    const double MAX_STEP = (pc.endparam() - pc.startparam())/20;
    const double dx_max = max( *max_element(dx, dx + num_internal_points),
  			      -*min_element(dx, dx + num_internal_points));
    //cout << dx_max << endl;
    const double fac = (dx_max < MAX_STEP) ? 1 : MAX_STEP/dx_max;
    
    // trying to update function using jacobian.  If this fails, use the fact
    // that the linear system is really a gradient, and employ steepest descent.
    const double orig_funval = get<0>(fun)(x, (int)num_internal_points);
    vector<double> x_tmp(x, x + num_internal_points);

    for (size_t i = 0; i != (size_t)num_internal_points; ++i) {
      x_tmp[i] += fac * dx[i];
      x_tmp[i] = max(x_tmp[i], pc.startparam());
      x_tmp[i] = min(x_tmp[i], pc.endparam());
    }
    for (size_t i = 0; i != (size_t)num_internal_points-1; ++i) {
      x_tmp[i] = min(x_tmp[i], x_tmp[i+1]);
    }
    double new_funval = get<0>(fun)(&x_tmp[0], (int)num_internal_points);
    if (new_funval < orig_funval) {
      copy(x_tmp.begin(), x_tmp.end(), x);
      return;
    }

    // if we got here, the use of the jacobian did not decrease our function
    // value.  Let us try conjugated gradient instead
    cout << "Resort to conjugated gradient" << endl;
    FunctorWrapper funwrap(fun, num_internal_points, pc.startparam(), pc.endparam());
    FunctionMinimizer<FunctorWrapper> f_minimizer((int)num_internal_points, funwrap, x, TOL);
    minimise_conjugated_gradient(f_minimizer);
    copy(f_minimizer.getPar(), f_minimizer.getPar() + num_internal_points, x);

  };


  // find optimal set of points
  return find_root(get<1>(fun), &p[1], num_internal_points, updater, TOL, ITER); 
}
  
};
// ============================================================================
// --------------------------------- OLD CODE ---------------------------------
// ============================================================================


  // UpdateFun updater = [&] (double* x, const double* const dx) {

  //   const double MAX_STEP = (pc.endparam() - pc.startparam())/20;
  //   const double dx_max = max( *max_element(dx, dx + num_internal_points),
  // 			      -*min_element(dx, dx + num_internal_points));
  //   //cout << dx_max << endl;
  //   const double fac = (dx_max < MAX_STEP) ? 1 : MAX_STEP/dx_max;
    
  //   // trying to update function using jacobian.  If this fails, use the fact
  //   // that the linear system is really a gradient, and employ steepest descent.
  //   const double orig_funval = get<0>(fun)(x, (int)num_internal_points);
  //   vector<double> x_tmp(x, x + num_internal_points);

  //   for (size_t i = 0; i != (size_t)num_internal_points; ++i) {
  //     x_tmp[i] += fac * dx[i];
  //     x_tmp[i] = max(x_tmp[i], pc.startparam());
  //     x_tmp[i] = min(x_tmp[i], pc.endparam());
  //   }
  //   for (size_t i = 0; i != (size_t)num_internal_points-1; ++i) {
  //     x_tmp[i] = min(x_tmp[i], x_tmp[i+1]);
  //   }
  //   double new_funval = get<0>(fun)(&x_tmp[0], (int)num_internal_points);
  //   if (new_funval < orig_funval) {
  //     copy(x_tmp.begin(), x_tmp.end(), x);
  //     return;
  //   }

  //   // if we got here, the use of the jacobian did not decrease our function
  //   // value.  Let us try a steepest descent instead.
  //   cout << "Resort to SD" << endl;
  //   const auto res = get<1>(fun)(x, num_internal_points);
  //   vector<double> dir(res.value);
  //   transform(dir.begin(), dir.end(), dir.begin(), [](double x) {return -x;});
  //   const double init_steplength = sqrt(norm2(dx, num_internal_points));
  //   const double max_steplength = compute_max_steplength(x, dx, num_internal_points,
  // 							 pc.startparam(), pc.endparam());

  //   cout << "Cur fval: " << get<0>(fun)(x, num_internal_points) << endl;
  //   // cout << "x is : ";
  //   // copy(x, x+num_internal_points, ostream_iterator<double>(cout, ", ")); cout << endl;
  //   x_tmp = directional_minimum(get<0>(fun),
  // 				x,
  // 				&dir[0],
  // 				num_internal_points,
  // 				init_steplength,
  // 				max_steplength,
  // 				TOL * (pc.endparam() - pc.startparam()));
  //   cout << "New fval: " << get<0>(fun)(&x_tmp[0], num_internal_points) << endl;
  //   // cout << "new x is : ";
  //   // copy(x_tmp.begin(), x_tmp.end(), ostream_iterator<double>(cout, ", ")); cout << endl;

    
  //   copy(x_tmp.begin(), x_tmp.end(), x);
  // };


  // ----------------------------------------------------------------------------
  // double compute_max_steplength(const double* const x,
  // 				const double* const dx,
  // 				const unsigned int dim,
  // 				const double startparam,
  // 				const double endparam)
  // // ----------------------------------------------------------------------------    
  // {
  //   vector<double> spacing(dim+1);
  //   for (unsigned int i = 0; i != dim+1; ++i) 
  //     spacing[i] =
  // 	(i==0)   ? x[0] - startparam :
  // 	(i==dim) ? endparam - x[dim-1] :
  // 	x[i] - x[i-1];

  //   vector<double> maxstep(dim);
  //   for (unsigned int i = 0; i != dim; ++i)
  //     maxstep[i] =
  // 	(dx[i] > 0) ? spacing[i+1] / dx[i] :
  // 	(dx[i] < 0) ? spacing[i] / abs(dx[i]) :
  // 	              numeric_limits<double>::infinity();
  //   return *min_element(maxstep.begin(), maxstep.end());
  // };

  // // ----------------------------------------------------------------------------
  // const pair<double, double> compute_bracket(const function<double(double)>& fun,
  // 					     const double init_guess,
  // 					     const double max_step)
  // // ----------------------------------------------------------------------------
  // {
  //   const double f0 = fun(0);

  //   double x2 = init_guess;
  //   double x1 = x2/2;
  //   double f2 = fun(x2);
  //   double f1 = fun(x1);

  //   // Since the function should have a negative derivative at origin, the following loop
  //   // SHOULD eventually exit.
  //   while (f1 > f0) {
  //     x2 = x1;
  //     f2 = f1;
  //     x1 = x2/2;
  //     f1 = fun(x1);
  //   }

  //   while (f2 < f1) {
  //     x2 = min(x2 * 2, max_step);
  //     f2 = fun(x2);

  //     if (x2==max_step) // we did not manage to bracket before reaching the boundary
  // 	break;
  //   }
  //   return {x1, x2};
  // };

  // // ----------------------------------------------------------------------------
  // const double golden_search(const function<double(double)>& fun,
  // 			     const pair<double, double>& bracket,
  // 			     const double tolerance)
  // // ----------------------------------------------------------------------------
  // {
  //   // @@ DEBUG PURPOSE
  //   const double f_start = fun(0);
  //   const double f_end = fun(bracket.second);
  //   // END DEBUG PURPOSE
  //   if (fun(bracket.second) < fun(bracket.first))
  //     return bracket.second; // for some reason, bracketing could not be achieved
    
  //   const double GOLD = 0.618034;
  //   const double GOLD_COMPL = 1 - GOLD;
  //   double x0 = 0;
  //   double x1, x2;
  //   double x3 = bracket.second;
  //   if (abs(x3 - bracket.first) > abs(bracket.first - x0)) {
  //     x1 = bracket.first;
  //     x2 = x1 + GOLD_COMPL * (x3 - x1);
  //   } else {
  //     x2 = bracket.first;
  //     x1 = x2 - GOLD_COMPL * (x2 - x0);
  //   }
  //   double f1 = fun(x1);
  //   double f2 = fun(x2);
  //   while (abs(x3-x0) > tolerance * (x1 + x2)) {
  //     if (f2 < f1) {
  // 	x0 = x1;
  // 	x1 = x2;
  // 	x2 = GOLD * x1 + GOLD_COMPL * x3;
  // 	f1 = f2;
  // 	f2 = fun(x2);
  //     } else {
  // 	x3 = x2;
  // 	x2 = x1;
  // 	x1 = GOLD * x2 + GOLD_COMPL * x0;
  // 	f2 = f1;
  // 	f1 = fun(x1);
  //     }
  //   }

  //   return (f1 < f2) ? x1 : x2;
  // };

  // // ----------------------------------------------------------------------------  
  // vector<double> directional_minimum(const RnToRFunction& fun,
  // 				     const double* const t,
  // 				     const double* const dir,
  // 				     const unsigned int dim,
  // 				     const double init_guess,
  // 				     const double max_length,
  // 				     const double tolerance)
  // // ----------------------------------------------------------------------------    
  // {
  //   // search for a minimum of function for parameter t + alpha * dir
  //   vector<double> tmp(dim, 0);
  //   function<double(double)> univariate_fun = [&](double alpha) {
  //     transform(t, t+dim, dir, tmp.begin(), [alpha] (double a, double b)->double {return a + alpha * b;});
  //     return fun(&tmp[0], (int)dim);
  //   };
    
  //   const pair<double, double> bracket = compute_bracket(univariate_fun, init_guess, max_length);
  //   const double alpha = golden_search(univariate_fun, bracket, tolerance);
  //   vector<double> res(dim);
  //   transform(t, t+dim, dir, res.begin(), [alpha](double a, double b) {return a + alpha * b;});
  //   return res;
  // }

  // // ----------------------------------------------------------------------------
  // double norm2(const double* const x, unsigned int num)
  // // ----------------------------------------------------------------------------
  // {
  //   double res = 0;
  //   for (unsigned int i = 0; i != num; ++i) {
  //     res = res + x[i] * x[i];
  //   }
  //   return res;
  // }
    
