#include "lu.h"
#include "GoTools/utils/LUDecomp.h"

using namespace std;
using namespace Go;

namespace {

// wrapper class used to convert a matrix stored column-wised in a contigous
// memory array into an object that the Go::LUDecomp function can deal with.
class Mat {
public:
  Mat(unsigned int N, double* data) : p_(N, data) {}
  struct Proxy {
    Proxy(unsigned int N, double* data) : n_(N), data_(data) {}
    double& operator[](int j) { return data_[j * n_ + i_]; }
    unsigned int n_;
    unsigned int i_;
    double* data_;
  };
  Proxy& operator[](int i) { p_.i_ = i; return p_;}
private:
  Proxy p_;
};

  
}; // anonymous namespace

namespace TesselateUtils {

// ============================================================================
// at the moment, this function is just a wrapper for the GoTools
// implementation, which requires a matrix object with the [][] operator.
bool lu(unsigned int N, double* const coefs, int* perm, bool& parity)
// ============================================================================  
{
  try {
    Mat m(N, coefs);
    LUDecomp(m, N, perm, parity);
    return true;
  } catch (exception& e) {
    // an exception thrown from LUDecomp means that system was singular
    return false;
  }
  
}
  
}; // end namespace TesselateUtils
