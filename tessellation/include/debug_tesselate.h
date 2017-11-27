#ifndef _DEBUG_TESSELATE_H
#define _DEBUG_TESSELATE_H

#include <vector>
#include <string>
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "clip_grid.h"
#include "common_defs.h"

using namespace TesselateUtils;
using namespace std;


void sample_polyhedron_energy(double xmin, double xmax, int num_x,
                              double ymin, double ymax, int num_y,
                              double zmin, double zmax, int num_z,
                              const Point3D* const bpoints,
                              const unsigned int num_bpoints,
                              const Triangle* const btris,
                              const unsigned int num_btris,
                              const double vdist,
                              string filename);

vector<int> identify_nondelaunay_tets(const vector<Point3D>& points,
                                      const vector<Tet>& tets,
                                      const uint num_bpoints);

vector<uint> interior_tets(const vector<Tet>& tets, const uint num_bpoints);
                               

#endif
