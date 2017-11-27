#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string> // stoi
#include <limits>

// #define SILENCE_EXTERNAL_WARNINGS 1
// #include "disable_warnings.h"
// #include <Eigen/Dense>
// #include <Eigen/Core>
// #include "reenable_warnings.h"

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/PointCloud.h"
#include "make_spacing_funs.h"
#include "find_root.h"

using namespace std;
using namespace Go;
//using namespace Eigen;

namespace {
  vector<double> define_parameters(double minparam, double maxparam, unsigned int num_intparams) {
    vector<double> result(num_intparams+2);
    result.front() = minparam;
    result.back()  = maxparam;
    for (size_t i = 1; i != num_intparams+1; ++i)
      result[i] = (maxparam-minparam)/(num_intparams + 1) * double(i);
    return result;
  }


  PointCloud<3> compute_pcloud(const SplineCurve& c, const vector<double>& pars)
  {
    vector<double> points(pars.size() * 3);
    for (size_t i = 0; i != pars.size(); ++i) {
      Point p;
      c.point(p, pars[i]);
      copy(p.begin(), p.end(), &points[3*i]);
    }
    //pcloud = PointCloud<3>(points.begin(), (int) points.size()/3);
    return PointCloud<3>(points.begin(), (int)points.size()/3);
  }
    
};

int main(int varnum, char* vararg[])
{
  const unsigned int num_int_points = (varnum > 1) ? std::stoi(vararg[1]) : 5;
  
  // make sample spline curve
  const SplineCurve c(4, 4, &(vector<double>({0, 0, 0, 0, 1, 1, 1, 1}))[0],
		      &(vector<double>({0,    0, 0,
			                0.25, 1, 0,
			                0.75, 0, 0,
			                1, 0, 0}))[0], 3);

  // define a set of parameters
  const vector<double> pars = define_parameters(c.startparam(),
						c.endparam(),
						num_int_points);
  
  const auto fun = make_curve_spacing_fun(c);

  // search for optimum
  UpdateFun updater = [&] (double* x, const double* const dx) {

    const double MAX_STEP = (c.endparam() - c.startparam())/20;
    const double dx_max = max( *max_element(dx, dx + num_int_points),
			      -*min_element(dx, dx + num_int_points));
    cout << dx_max << endl;
    const double fac = (dx_max < MAX_STEP) ? 1 : MAX_STEP/dx_max;
    
    for (size_t i = 0; i != num_int_points; ++i) {
      x[i] += fac * dx[i];
      x[i] = max(x[i], c.startparam());
      x[i] = min(x[i], c.endparam());
    }

    // for (size_t i = 1; i != num_int_points; ++i) {
    //   x[i] = max(x[i], x[i-1]);
    // }
    for (size_t i = 0; i != num_int_points-1; ++i) {
      x[i] = min(x[i], x[i+1]);
    }
  };

  const vector<double> new_pars = find_root(get<1>(fun), &pars[1], num_int_points, updater, 1e-9, 20);

  // Text reporting
  cout << "Original vector: \n";
  copy(pars.begin()+1, pars.end()-1, ostream_iterator<double>(cout, ", "));
  cout << endl;
  cout << "Associated function value: " << get<0>(fun)(&pars[1], num_int_points) << endl;
  
  cout << "New vector: \n" ;
  copy(new_pars.begin(), new_pars.end(), ostream_iterator<double>(cout, ", "));
  cout << endl;
  cout << "Associated function value: " << get<0>(fun)(&new_pars[0], num_int_points) << endl;
  cout << endl;
  // Save results for plotting

  ofstream os("result.g2");
  c.writeStandardHeader(os);
  c.write(os);

  const PointCloud<3> pcloud = compute_pcloud(c, pars);
  pcloud.writeStandardHeader(os);
  pcloud.write(os);

  const PointCloud<3> pcloud_new = compute_pcloud(c, new_pars);
  pcloud_new.writeStandardHeader(os);
  pcloud_new.write(os);
  
};
