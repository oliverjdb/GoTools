#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include <iterator>

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/PointCloud.h"

#include "distributeND.h"

using namespace std;
using namespace Go;

namespace {
  void test_straight_line(int varnum, char* vararg[]);
  void test_spline_curve(int varnum, char* vararg[]);

  PointCloud<3> compute_pcloud(const SplineCurve& c, const vector<double>& pars);

  double parcurve_distfun(const double* p1, const double* p2, double* grad,
			  double* jac, const SplineCurve& c);
  
};

// ----------------------------------------------------------------------------
int main(int varnum, char* vararg[])
// ----------------------------------------------------------------------------
{
  //test_straight_line(varnum, vararg);
  test_spline_curve(varnum, vararg);
};

namespace {

// ----------------------------------------------------------------------------  
void test_spline_curve(int varnum, char* vararg[])
// ----------------------------------------------------------------------------
{
  //const int num_pts = stoi(vararg[1]);
  
  // make sample spline curve
  const SplineCurve c(4, 4, &(vector<double>({0, 0, 0, 0, 1, 1, 1, 1}))[0],
		      &(vector<double>({0,    0, 0,
			                0.25, 1, 0,
			                0.75, 0, 2,
			                1, 0, 0}))[0], 3);

  //const vector<double> input_points {0.1, 0.2, 0.3, 0.11, 0.111, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.9, 0.91, 0.92, 0.93, 0.933};
  //const vector<double> input_points {0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.50, 0.55, 0.60, 0.66, 0.7, 0.75, 0.80, 0.85, 0.90, 0.93, 0.95, 0.15};

  vector<double> input_points; 
  for (int i = 0; i != 100; ++i) {
    double tmp = 0.5 + i*0.001;
    input_points.push_back(tmp*tmp*tmp);
  }
  
  vector<double> result =
    distribute_points_1D(&input_points[0], (int)input_points.size(),
			 {c.startparam(), c.endparam()},
			 [&c] (const double* p1,
			       const double* p2,
			       double* grad,
			       double* jac) {
			   return parcurve_distfun(p1, p2, grad, jac, c);
			 });			 
  cout << "Result: " << endl;
  sort(result.begin(), result.end());
  copy(result.begin(), result.end(), ostream_iterator<double>(cout, ", "));
  cout << endl << endl;

  // saving
  ofstream os("result.g2");
  c.writeStandardHeader(os);
  c.write(os);

  const PointCloud<3> pcloud = compute_pcloud(c, result);
  pcloud.writeStandardHeader(os);
  pcloud.write(os);
  
} 
  
  
// ----------------------------------------------------------------------------
void test_straight_line(int varnum, char* vararg[])
// ----------------------------------------------------------------------------
{
  int num_pts = varnum - 3;
  const double min_bnd = stod(vararg[1]);
  const double max_bnd = stod(vararg[2]);
  
  vector<double> points(num_pts);
  for (int i = 0; i != num_pts; ++i)
    points[i] = stod(vararg[i+3]);
		     
  vector<double> result =
    distribute_points_1D(&points[0], num_pts, {min_bnd, max_bnd},
			 [] (const double* p1, const double* p2, double* der, double* dder) {
			   if (der) {
			     der[0] = (p1 > p2) ?  1 : -1;
			     der[1] = (p1 > p2) ? -1 :  1;
			   }
			   if (dder) {
			     dder[0] = 0;
			     dder[1] = 0;
			     dder[2] = 0;
			   }
			   return *p2 - *p1;
			 });

  cout << "Result: " << endl;
  copy(result.begin(), result.end(), ostream_iterator<double>(cout, ", "));
  cout << endl << endl;
}

// ----------------------------------------------------------------------------
PointCloud<3> compute_pcloud(const SplineCurve& c, const vector<double>& pars)
// ----------------------------------------------------------------------------
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
  
// ----------------------------------------------------------------------------
double parcurve_distfun(const double* p1, const double* p2, double* grad,
			double* jac, const SplineCurve& c)
// ----------------------------------------------------------------------------
{
  vector<Point> pts1(3), pts2(3);
  c.point(pts1, *p1, 2);
  c.point(pts2, *p2, 2);

  const Point d = pts2[0] - pts1[0];
  const double dist = d.length() > 0 ? d.length() : 1;
  const double dist3 = dist * dist * dist;
  const double d_dp1 = d * pts1[1];
  const double d_dp2 = d * pts2[1];
  const double d_ddp1 = d * pts1[2];
  const double d_ddp2 = d * pts2[2];

  if (grad) {
    grad[0] = - d_dp1 / dist;
    grad[1] =   d_dp2 / dist;
  }
  if (jac) {
    jac[0] = (pts1[1].length2() - d_ddp1) / dist - (d_dp1 * d_dp1 / dist3);// d2/du2
    jac[1] = (pts2[1].length2() + d_ddp2) / dist - (d_dp2 * d_dp2 / dist3); // d2/dv2
    jac[2] = -(pts1[1] * pts2[1]) / dist         + (d_dp1 * d_dp2 / dist3); // d2/dudv
  }
  
  return dist;
}
 
}; // end anonymous namespace
