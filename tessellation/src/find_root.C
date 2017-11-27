#include <stdexcept>
#include <iostream>

#define SILENCE_EXTERNAL_WARNINGS 1
#include "disable_warnings.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include "reenable_warnings.h"

#include "make_spacing_funs.h"

using namespace std;
using namespace Go;
using namespace Eigen;

namespace Go {

std::vector<double>
find_root(const RnToRnFunction& fun,
	  const double* const start,
	  unsigned int dim,
	  const UpdateFun& ufun,
	  double tol,
	  unsigned int max_iter)
{
  MatrixXd jacobian(dim, dim);
  VectorXd d(dim), dv(dim), v(dim);
  copy(start, start+dim, v.data());

  for (size_t i = 0; i != max_iter; ++i) {
    cout << i << endl;
    // compute new values and jacobian
    const ValAndJac val = fun(v.data(), dim);
    copy(val.value.begin(), val.value.end(), d.data());
    copy(val.jacobian.begin(), val.jacobian.end(), jacobian.data());

    // check if we are at a root.  If not, do an update
    cout << "Error: " << d.cwiseAbs().maxCoeff() << endl;
    if (d.cwiseAbs().maxCoeff() < tol) {
      cout << "Converged in " << i << " iterations." << endl;
      break;
    }
    
    // compute update, according to Newton
    dv = -1 * jacobian.colPivHouseholderQr().solve(d);

    // update search point according to provided policy
    ufun(v.data(), dv.data());
  }
  
  if (d.cwiseAbs().maxCoeff() > tol)
    throw runtime_error("Failed to converge in 'find_root'");

  return vector<double> {v.data(), v.data() + dim};
}
  
  
}; // end namespace Go
