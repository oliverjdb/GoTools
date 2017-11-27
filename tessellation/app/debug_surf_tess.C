#include <fstream>
#include <iostream>
#include <vector>
#include <iterator>
#include "tesselate_utils.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/CPUclock.h"
#include "interpoint_distances.h"
#include "polyhedral_energies.h"
#include "triangulate_domain.h"
#include "SimplePolyhedronTesselation.h"
#include "GoParametricTesselableVolume.h"
#include "tesselate_polyhedron.h"
#include "fit_points_to_plane.h"
#include "basic_intersections.h"
#include "clip_grid.h"

#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"
//#include "GoTools/compositemodel/CompositeModelFileHandler.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "tesselate_parametric_volume.h"

using namespace Go;
using namespace std;
using namespace TesselateUtils;

typedef shared_ptr<const ParamSurface> SurfPtr;

int main() {


  const vector<double> knots = {0, 0, 0, 0, 1, 1, 1, 1};
  vector<double> coefs;
  for(int j = 0; j != 4; ++j)
    for (int i = 0; i != 4; ++i) {
      coefs.push_back(i/3.0);
      coefs.push_back(j/3.0);
      coefs.push_back(0.0);
    }
          
  SurfPtr s(new SplineSurface(4, 4, 4, 4, &knots[0], &knots[0], &coefs[0], 3));

  //vector<Point2D> bnd { {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0} };
  vector<Point2D> bnd { {0.5, 0.0}, {1.0, 0.5}, {0.5, 1.0}, {0.0, 0.5} };  
  
  Mesh2D m = tesselateParametricSurface(s, &bnd[0], 4, 0.1);
  
  
};
