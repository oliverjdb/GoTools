#ifndef _LU_H
#define _LU_H

namespace TesselateUtils {

// ============================================================================
// compute LU decomposition of an N x N matrix, where column-wise coefficients
// are stored in 'coefs' (will contain the decomposition upon return).  The
// row-wise permutation is returned in 'perm', which must point to an N-sized
// array.  'paritty' reports whether the number of row interchanges were pair or
// impair.  Return value indicates whether decomposition was successful or not.
bool lu(unsigned int N, double* const coefs, int* perm, bool& parity);
// ============================================================================
  
}; // end namespace TesselateUtils

#endif
