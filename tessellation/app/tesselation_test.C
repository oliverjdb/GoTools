#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm>

#define SILENCE_EXTERNAL_WARNINGS 1
#include "disable_warnings.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include "reenable_warnings.h"

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/PointCloud.h"

using namespace std;
using namespace Go;
using namespace Eigen;

struct Res {
  double val;
  vector<double> d1;
  vector<double> d2_diag;
  vector<double> d2_subdiag;
};

vector<double> prepostpend(double first, double last,
			   const double* const mid_start, size_t mid_size)
{
  vector<double> result(mid_size+2);
  result.front() = first;
  result.back() = last;
  copy(mid_start, mid_start + mid_size, &result[1]);
  return result;
}
			   

// The values in u should be monotonously increasing, and be strictly between
// c.startparam() and c.endparam().
Res splinePtFun(const SplineCurve& c, const double* const u, size_t num)
{
  // making full parameter vector, including start and end parameters
  const vector<double> u_full = prepostpend(c.startparam(), c.endparam(), u, num);

  // vector for storing points and derivatives
  vector<vector<Point>> points(u_full.size(), vector<Point>(3));

  // evaluate all points and their first and second derivatives
  auto to = points.begin();
  for (auto from = u_full.begin(); from != u_full.end(); ++from, ++to)
    c.point(*to, *from, 2);

  // compute squared distances
  double val = 0;
  for (auto p1 = points.begin(), p2 = p1+1; p2 != points.end(); ++p1, ++p2) {
    val = val + (*p1)[0].dist2((*p2)[0]);
  }

  // computing first derivative
  vector<double> d1(num);
  for (size_t i = 1; i != num+1; ++i) {
    d1[i-1] = (4 * points[i][0] - 2 * (points[i+1][0] + points[i-1][0])) * points[i][1];
  }

  // computing second derivative
  vector<double> d2_diag(num), d2_subdiag(num-1);
  for (size_t i = 1; i != num+1; ++i) {
    d2_diag[i-1] =
      (4 * points[i][0] - 2 * (points[i+1][0] + points[i-1][0])) * points[i][2] +
      4 * points[i][1].length2();
    if (i>1) d2_subdiag[i-2]   = -2 * points[i-1][1] * points[i][1];
  }
  return {val, d1, d2_diag, d2_subdiag};
};

Matrix<double, Dynamic, Dynamic> construct_jacobian(const Res& res) {
  const size_t dim = res.d1.size();
  Matrix<double, Dynamic, Dynamic> m(dim, dim);
  m *= 0;
  for (size_t i = 0; i != (size_t)dim; ++i) {
    m(i, i) = res.d2_diag[i];
    if (i>0)     m(i, i-1) = res.d2_subdiag[i-1];
    if (i<(size_t)dim-1) m(i, i+1) = res.d2_subdiag[i];
  }
  return m;
}

// ----------------------------------------------------------------------------

int main() {

  SplineCurve c(4, 4, &(vector<double>({0, 0, 0, 0, 1, 1, 1, 1}))[0],
		&(vector<double>({0,    0, 0,
			0.25, 1, 0,
			0.75, 0, 0,
			1, 0, 0}))[0], 3);

  ofstream os("result.g2");

  //vector<double> pars {0, 0.15, 0.5, 0.75, 1};
  //vector<double> pars {0, 0.25, 0.5, 0.75, 1};
  //vector<double> pars {0, 0.2, 0.4, 0.6, 0.8, 1};
  vector<double> pars {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

  vector<double> points(pars.size()*3);
  for (size_t i = 0; i != pars.size(); ++i) {
    Point p;
    c.point(p, pars[i]);
    copy(p.begin(), p.end(), &points[3*i]);
  }
  
  PointCloud<3> pcloud(points.begin(), (int)points.size()/3);

  // compute cost function
  auto res = splinePtFun(c, &pars[1], pars.size()-2);
  cout << "Cost function is: " << res.val << endl;
  cout << "Derivative is: \n";
  copy(res.d1.begin(), res.d1.end(), ostream_iterator<double>(cout, ", "));
  cout << endl << endl;
  cout << "Second derivative:\n";
  copy(res.d2_diag.begin(), res.d2_diag.end(), ostream_iterator<double>(cout, ", "));
  cout << endl;
  copy(res.d2_subdiag.begin(), res.d2_subdiag.end(), ostream_iterator<double>(cout, ", "));
  cout << endl;

  // making matrix of second derivative
  Matrix<double, Dynamic, Dynamic> m = construct_jacobian(res);

  // const int dim = (int)pars.size()-2;
  // Matrix<double, Dynamic, Dynamic> m(dim, dim);
  // m *= 0;
  // for (size_t i = 0; i != (size_t)dim; ++i) {
  //   m(i, i) = res.d2_diag[i];
  //   if (i>0)     m(i, i-1) = res.d2_subdiag[i-1];
  //   if (i<(size_t)dim-1) m(i, i+1) = res.d2_subdiag[i];
  // }

  // vector of parameter points
  size_t dim = res.d1.size();
  VectorXd v(dim); for (size_t i = 0; i != (size_t)dim; ++i) v(i) = pars[i+1];

  // Vector of derivative
  VectorXd d(dim); for (size_t i = 0; i != (size_t)dim; ++i) d(i) = res.d1[i];

  VectorXd dv;
  const double MAX_DIFF = 0.05;

  for (int iter = 1; iter != 8; ++iter) {

    cout << "====iteration " << iter << "====\n";
    
    // Computing update
    dv = -1 * m.colPivHouseholderQr().solve(d);

    // rescaling update to ensure no change is larger than MAX_DIFF
    double max_elem = abs(dv(0));
    for (size_t i = 1; i != dim; ++i)
      max_elem = max(max_elem, abs(dv(i)));
    if (max_elem > MAX_DIFF) {
      double fac = MAX_DIFF / max_elem;
      dv *= fac;
    }
    cout << "Matrix is: \n" << m << endl;

    cout << "Update is: \n" << dv << endl;

    VectorXd v_new = v+dv;
    // bounds checks
    v_new(0) = std::max(v_new(0), double(0));
    v_new[dim-1] = std::min(v_new[dim-1], double(1));
    for (size_t i = 1; i != dim; ++i)
      v_new[i] = max(v_new[i], v[i-1]);
    for (size_t i = 0; i != dim-1; ++i)
      v_new[i] = min(v_new[i], v[i+1]);

    cout << "New parameters are: \n" << v_new << endl;

    // computing new gradient
    Res res_new = splinePtFun(c, v_new.data(), dim);
    VectorXd d_new(dim);
    copy(res_new.d1.begin(), res_new.d1.end(), d_new.data());
    cout << "New gradient is: \n" << d_new << endl;

    v = v_new;
    d = d_new;
    m = construct_jacobian(res_new);
  }


  vector<double> points_new((size_t)dim * 3);
  for (size_t i = 0; i != (size_t)dim; ++i) {
    Point p;
    c.point(p, v(i));
    copy(p.begin(), p.end(), &points_new[3*i]);
  }
  
  PointCloud<3> pcloud_new(points_new.begin(), (int)points_new.size()/3);



  
  //m.diagonal();
  //  DiagonalMatrix<double, Dynamic> m(&(res.d2_diag[0]), 3);
  
  c.writeStandardHeader(os);
  c.write(os);
  pcloud.writeStandardHeader(os);
  pcloud.write(os);
  pcloud_new.writeStandardHeader(os);
  pcloud_new.write(os);

  os.close();

  return 0;
  
}
