#include <algorithm>
#include <iostream>
#include <cstdlib>
#include "tesselate_polyhedron.h"

using namespace std;
using namespace TesselateUtils;

// ----------------------------------------------------------------------------
int main(int varnum, char** vararg)
// first argument: target distance (double)
// second argument: number of polygon corners (int)
// following arguments: polygon corner coordinates
// ----------------------------------------------------------------------------
{
  if (varnum <= 2) {
    cout << "You need to supply arguments. \n";
    cout << "Try the following argument set, for instance:\n";
    cout << "app/polyhedron_tesselation_test 0.05 5 0 0 1 0 1 2 0.5 1 0 1" << endl;
    return 0;
  }

  const double l = atof(vararg[1]);
  const int num_corners = atoi(vararg[2]);

  // read corners
  vector<Point2D> corners;
  for (int i = 0; i != num_corners; ++i) 
    corners.emplace_back(Point2D {atof(vararg[2*i + 3]), atof(vararg[2*i+4])});

  const auto result = tesselatePolygon2D(&corners[0], (unsigned int)corners.size(), l);

  cout << "Points: " << endl;
  for (const auto p : result.points)
    cout << p[0] << " " << p[1] << '\n';
  cout << endl << "Tris: " << endl;
  for (const auto t : result.tris)
    cout << t[0] << " " << t[1] << " " << t[2] << '\n';
  
  
  // for(auto p = result.begin(); p != result.end(); ++p)
  //   cout << (*p)[0] << " " << (*p)[1] << '\n';
  
  return 0;
};

  
