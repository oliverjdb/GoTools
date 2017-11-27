#include <cmath>
#include "ParametricObjectEnergyFunctor.h"
#include "GoTools/utils/GeneralFunctionMinimizer.h"
#include "tesselate_parametric_volume.h"

using namespace TesselateUtils;
using namespace std;
using namespace Go;

namespace {

vector<Point> evaluate_points(const shared_ptr<const ParamCurve> pc,
                              const vector<double>& t);
vector<double> compute_distances(const vector<Point> points);
vector<double> define_parvec(double minparam,
                             double maxparam,
                             unsigned int num_intparams);
vector<double> adjust_parameters_to_geometry(const shared_ptr<const ParamCurve> pc,
                                             const vector<double>& t);
void optimize_interior_points(const shared_ptr<const ParamCurve> pc,
                              double radius,
                              double* const par,
                              uint num_par);
                              
  
}; // end anonymous namespace 

namespace TesselateUtils {

// ----------------------------------------------------------------------------
vector<double> tesselateParametricCurve(const shared_ptr<const ParamCurve> pc,
                                        const double vdist)
// ----------------------------------------------------------------------------
{
  // estimate lenght of curve, in order to determine the necessary number of
  // inner points
  const double DIST_FAC = 1.5; //1.5; // properly adjust radius of energy function
  const int NUM_SAMPLES = 20;
  const int N = max(0,
                    (int)floor(pc->estimatedCurveLength(NUM_SAMPLES) / vdist) - 1);

  // initial guess of parameters (we try to space them as evenly as possible)
  vector<double> par = adjust_parameters_to_geometry(pc,
                                                     define_parvec(pc->startparam(),
                                                                   pc->endparam(), N));
  if (!par.empty())
    optimize_interior_points(pc, vdist * DIST_FAC, &par[1], (uint)par.size()-2);  
  
  return par; 
}
  

// // ----------------------------------------------------------------------------
// double estimateCurveLength(const shared_ptr<const ParamCurve> pc,
//                            const unsigned int num_samples)
// // ----------------------------------------------------------------------------  
// {
//   const vector<double> par =
//     adjust_parameters_to_geometry(pc, define_parvec(pc->startparam(),
// 						    pc->endparam(),
// 						    num_samples-2));
//   Point p, pnew;
//   double accum = 0;
//   pc->point(p, par[0]);
//   for (size_t i = 1; i != par.size(); ++i) {
//     pc->point(pnew, par[i]);
//     accum += pnew.dist(p);
//     p.swap(pnew);
//   }
//   return accum;

// }

  
}; // end namespace TesselateUtils

namespace {

// ----------------------------------------------------------------------------  
vector<double> adjust_parameters_to_geometry(const shared_ptr<const ParamCurve> pc,
                                             const vector<double>& t)
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

// ----------------------------------------------------------------------------
vector<Point> evaluate_points(const shared_ptr<const ParamCurve> pc,
                              const vector<double>& t)
// ----------------------------------------------------------------------------
{
  vector<Point> result(t.size(), Point(3));
  for (size_t i = 0; i != t.size(); ++i)
    pc->point(result[i], t[i]);
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
vector<double> define_parvec(double minparam,
                             double maxparam,
                             unsigned int num_intparams)
// ----------------------------------------------------------------------------
{
  vector<double> result(num_intparams+2);
  result.front() = minparam;
  result.back()  = maxparam;
  for (size_t i = 1; i != num_intparams+1; ++i)
    result[i] = minparam + (maxparam-minparam)/(num_intparams + 1) * double(i);
  return result;
}

// ----------------------------------------------------------------------------  
void optimize_interior_points(const shared_ptr<const ParamCurve> pc,
                              double radius,
                              double* const par,
                              uint num_par)
// ----------------------------------------------------------------------------
{
  // setting up function to minimize
  const array<Point1D, 2> boundary = {pc->startparam(), pc->endparam()};
  auto efun = ParamCurveEnergyFunctor(pc, {&boundary[0]}, radius, num_par);

  // setting up function minimizer
  Go::FunctionMinimizer<ParamCurveEnergyFunctor>
    funcmin(num_par, efun, par, 1e-6);//1e-1); // @@ to lax tolerance?

  const double STOPTOL = 1e-6; // @@

  // Do the minimization
  Go::minimise_conjugated_gradient(funcmin, STOPTOL);

  // copy results back
  copy(funcmin.getPar(), funcmin.getPar() + num_par, par);
  
}
  
  
}; // end anonymous namespace 

