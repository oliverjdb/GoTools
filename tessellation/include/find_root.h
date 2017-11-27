#ifndef _FINDROOT_H
#define _FINDROOT_H

#include <vector>
#include <functional>

namespace Go {

// ----------------------------------------------------------------------------

// Class representing the value and jacobian of a multi-parameter, multi-valued function
struct ValAndJac
{
  std::vector<double> value; // of size N
  std::vector<double> jacobian; // of size NxN (column-wise storage)
};
  
// ----------------------------------------------------------------------------

// Function from Rn to Rn, returning a ValAndJac.  First argument is a pointer
// to the function arguments, second argument is the number of function arguments
typedef std::function<ValAndJac(const double* const, unsigned int)> RnToRnFunction;

// ----------------------------------------------------------------------------

// Function used internally to update search point
// First argument: old point (values will be updated)
// Second argument: proposed update
// Third argument: number of components (dimension of point)
typedef std::function<void(double*, const double* const)> UpdateFun;
  
// Function searching for root
std::vector<double>
find_root(const RnToRnFunction& fun,
	  const double* const start,
	  unsigned int dim,
	  const UpdateFun& ufun,
	  double tol = 1e-9,
	  unsigned int max_iter = 12);

  
}; //end namespace Go

#endif // _FINDROOT_H
