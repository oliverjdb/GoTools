#include "debug_tesselate.h"
#include "tesselate_utils.h"
#include "polyhedral_energies.h"
#include <fstream>

using namespace TesselateUtils;
using namespace std;

// ----------------------------------------------------------------------------
void sample_polyhedron_energy(double xmin, double xmax, int num_x,
                              double ymin, double ymax, int num_y,
                              double zmin, double zmax, int num_z,
                              const Point3D* const bpoints,
                              const unsigned int num_bpoints,
                              const Triangle* const btris,
                              const unsigned int num_btris,
                              const double vdist,
                              string filename)
// ----------------------------------------------------------------------------
{
  vector<Point3D> sample_points = generate_grid_3D(Point3D {xmin, ymin, zmin},
                                                   Point3D {xmax, ymax, zmax},
                                                   (uint) num_x,
                                                   (uint) num_y,
                                                   (uint) num_z,
                                                   false);

  // computing values
  vector<ValAndDer<Point3D>> vals;
  transform(sample_points.begin(), sample_points.end(), back_inserter(vals),
            [&] (const Point3D& p) { return polyhedron_energy(bpoints, num_bpoints,
                                                              btris, num_btris,
                                                              &p, 1, vdist); });
  
  // saving results
  ofstream os(filename.c_str());

  os << "A - samples: " << endl;
  for (auto p : sample_points)
    os << p[0] << " " << p[1] << " " << p[2] << '\n';
  
  os << '\n' << "B - values: " << endl;
  for (auto v : vals)
    os << v.val << ' ',

  os << '\n' << "C - derivatives: " << endl;
  for (auto v : vals) {
    for (auto d : v.der)
      os << d[0] << ' ' << d[1] << ' ' << d[2] << '\n';
    os << "\n\n";
  }
  
  os.close();
}

// ----------------------------------------------------------------------------
int empty_interior(const Tet& t, const vector<Point3D>& points, const uint num_bpoints)
// ----------------------------------------------------------------------------
{
  Point3D c; // center
  double r2; // radius squared
  fitting_sphere(points[t[0]], points[t[1]], points[t[2]], points[t[3]], c, r2);

  //for (uint i = 0; i != points.size(); ++i)
  for (uint i = num_bpoints; i != points.size(); ++i)  
    if ((i != t[0]) && (i != t[1]) && (i != t[2]) && (i != t[3]) &&
        (dist2(points[i], c) < r2))
      return i;
  return -1;
}

// ----------------------------------------------------------------------------
vector<int> identify_nondelaunay_tets(const vector<Point3D>& points,
                                      const vector<Tet>& tets,
                                      const uint num_bpoints)
// ----------------------------------------------------------------------------
{
  vector<int> result(tets.size());
  transform(tets.begin(), tets.end(), result.begin(),
            [&] (const Tet& tet) { return empty_interior(tet, points, num_bpoints);});

  return result;
}
  
// ----------------------------------------------------------------------------
vector<uint> interior_tets(const vector<Tet>& tets, const uint num_bpoints)
// ----------------------------------------------------------------------------
{
  vector<uint> result(tets.size());
  for (uint i = 0; i != (uint)tets.size(); ++i)
    result[i] = ( (tets[i][0] >= num_bpoints) &&
                  (tets[i][1] >= num_bpoints) &&
                  (tets[i][2] >= num_bpoints) &&
                  (tets[i][3] >= num_bpoints) );
  return result;
}
