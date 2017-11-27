#ifndef TESSELATE_DEBUG_
#define TESSELATE_DEBUG_

#include <string>
#include <fstream>

#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ParamCurve.h"

using namespace Go;
using namespace std;

// ----------------------------------------------------------------------------  
// @@ Debug purposes
void store_points_and_curve(const ParamCurve& pc,
			    double* t,
			    unsigned int num_points,
			    const string& filename);

#endif
